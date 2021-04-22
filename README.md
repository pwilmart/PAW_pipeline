# Update 6/30/2020
Made some changes to `make_PAW_TXT_from_PD2.x.py` to support Mascot search exports and some output changes in newer PD 2.x exports. There are  also some output log information changes for some pipeline steps. Added support for real time search SPS-MS3 raw files and 16-plex TMT reagents to `msconvert_GUI.py`. Added support for 16-plex TMT reagents and high res MS2 fragment ions to `comet_GUI.py`.

## Update 8/22/2019
Prior version of `sqt_converter.py` had a serious bug that resulted in a large amount of data to be dropped. Please discard any scripts downloaded between 8/6/2019 and 8/22/2019 and replace with new versions.

## Update 8/6/2019
New version of `sqt_converter.py` greatly speeds up processing of semi-trpytic searches. There are some changes to `PAW_lib.py` that are also needed. New versions of these two scripts are now in the repository. The conversion step creates peptide lookup dictionaries by doing an *in silico* digestion of the FASTA file. This can take a few minutes upfront, depending on database size. The dictionaries can be pretty large.

## Update 5/27/2019

Did some reformatting of summary tables (`results_files` folder) to make them easier to read in R scripts. The new table names will end in `_9` to help keep versions clear. It is easy to skip lines at the top of files when reading data into R, so some meta data (if any) will be in the first 4 rows. The table header line will be in row 5 and there will be nothing after the main table. The protein grouping script used to put a few things (coverage, length, MW, description) in parentheses because these may not have all been the same for all proteins that were grouped. That created columns with mixtures of numerical and text data. R generally does not like character/text data, so values from the primary protein are now used. There were also minor changes to some log files to collect additional meta data.

A new `PAW_table_descriptions_9.txt` has been added to explain the columns in all of the tables written by the PAW pipeline. Most of those are results tables, but there are also top hit summaries and some TMT data tables. Depending on pipeline options and data, not all tables will be written for each analysis. The description file lists all potential tables.

---

## List of scripts

### Comet-based pipeline scripts (in general order):

- **msconvert_GUI.py** - Converts Thermo RAW files into MS2 format files
- **comet_GUI.py** - Configures Comet params files
- **histo_GUI.py** - Interactive score threshold setting tool
- **PAW_results.py** - parsimonious protein inference and grouping

### If the data was a TMT experiment:

- **add_TMT_intensities.py** - aggregates PSM reporter ions into protein totals
- [optional] **pandas_TMT_IRS_norm.py** - performs IRS normalization for multi-TMT experiments

### If processing Proteome Discoverer exports:

- **make_PAW_TXT_from_PD1.4.py** - makes PAW files from PD 1.4 exports
- **make_PAW_TXT_from_PD2.x.py** - makes PAW files from PD 2.x exports
- **PD1.4_TMT_phospho_processer.py** - phosphopeptide TMT summarization

> **See "PAW_pipeline_for_PD_exports.pptx" for help with processing PD exports.**

### Support scripts:

- **PAW_lib.py** - main library of classes and functions
- **PAW_protein_grouper.py** - protein grouping script (called by PAW_results.py)

# PAW_pipeline

A Comet-based, best practices proteomics pipeline.

* [Introduction](#introduction)
* [Background](#background)
* [Best practices](#best-practices)
* [Pipeline Steps](#pipeline-steps)
* [Dependancies](#dependancies)
* [References](#references)

## Introduction
The PAW pipeline [1] consists of open source Python scripts to process Thermo mass spectrometer RAW files in a best-practices shotgun proteomics analysis pipeline. The open source ProteoWizard toolkit [2, 3] is used to convert scans in the RAW files into more easily readable formats. High-resolution isobaric labeling data (10- or 11-plex TMT reagents) can also be handled. Scripts simplify configuring and running MSConvert, convert the MS2 scans into appropriate file formats for database search, and extract the TMT reporter ion peak heights.  

The open source Comet database search engine [4] is used to identify peptides. Scripts simplify setting relevant search parameters, launching the Comet searches, and creating top-hit summary files. Other scripts [5] assist in downloading high-quality, up-to-date protein FASTA databases and adding decoy sequences and common contaminants.

Comet results are accurately filtered to desired false discovery rates (FDRs) using the target decoy method [6, 7] and proper use of accurate mass [8]. A GUI script generates accurate-mass-conditioned score histograms and FDR tables. The data are subclassed by peptide charge state and modification state. The scores are a linear discriminant function transformation similar to PeptideProphet [9] and Scaffold [10]. Score thresholds are set interactively in each data subclass. Inspection of score histograms provides integrated quality control.

Protein inference follows basic parsimony rules [11]. Multiple fractions per biological sample and multiple biological samples per experiment are supported. Protein inference is performed experiment-wide to provide robust results. Optional extended parsimony logic [12] is also implemented to group highly homologous proteins and maximize the quantitative information content of peptides. Reporter ion quantitative information in TMT experiments are aggregated to protein expression totals. Final protein and peptide reports are complete, carefully curated tab-delimited text files.

## Background
Analyzing bottom-up (a.k.a. shotgun) proteomics data can be a complicated, lengthy ordeal. The processing pipelines have many steps and each step can have several choices. There are few practical guidelines for making better choices to navigate this minefield. This results in far too great a diversity of pipelines in use, and most will have less-than-ideal performance.

This state of affairs is illustrated by one of the first steps of selecting a FASTA protein database file for use by the search engine. This is a more complicated topic than might be guessed [12]. There are often several protein database choices for many commonly studied organisms (and even some genomes that have not yet been sequenced). The protein database complexity and apparent size can vary considerably among the choices, although larger is probably not better [13]. Even integrated database sources like [UniProt](http://www.uniprot.org/) offer several versions of protein databases (Swiss-Prot, TrEMBL, canonical, canonical with isoforms, etc.), and selecting among them requires some expertise and use of [dedicated tools](https://github.com/pwilmart/fasta_utilities.git) [5].  

Picking a search engine program to identify likely peptide sequences associated with tandem mass spectra, known as peptide spectrum matches (PSMs), can be even more challenging. There are commercial products like Mascot, Proteome Discoverer, and Byonic that can be quite expensive. Alternatively, there are widely used freely available options such as X!Tandem and MaxQuant. There are also many less commonly used open source options, and opinions run high on the relative merits of the various search engines. Comet [4], written and maintained by Jimmy Eng, is an open source version of SEQUEST [14]. Comet is free, fast, sensitive, and a little simpler to configure than other options; this is what the PAW pipeline [1] uses. PAW are my initials, but you can make up something like Proteomic Analysis Workflow if you prefer.

The basic method of identifying peptides from tandem mass spectra (MS2 scans) has not really changed too much in 25 years [14]. It basically involves noise filtering and normalizing an MS2 scan, then comparing that to mass-filtered candidate peptide spectra from a theoretical digest of a protein database. Despite this clever leveraging of genomic sequencing (the protein sequences), the challenge has always been in deciding if a particular PSM is correct or not. Early heuristic approaches [15] led to basic classifiers [9], and even machine learning methods [16]. A big step forward came when decoy databases were used to eavesdrop on random noise scores, as popularized by the Gygi lab [6, 7]. More recent instrumental advances have added high-resolution and accurate mass to the equation [8].

With more robust statistical methods for controlling PSM errors, confident lists of identified **peptides** present in digests of complex protein mixtures can be reliably determined. The next issue is to determine a confident list of the **proteins** present in the samples from the identified peptides. This is known as the protein inference problem [11] and it persists despite significant advances in genomic sequencing (maybe that has made the problem even worse). Guidelines for parsimonious protein identifications [17] have now been widely used for many years. The basic parsimony rules are outlined in [11], but may need to be extended [12] for large datasets now routinely acquired.

## Best Practices
As proteomics has matured, there have been many analysis ideas that have come and gone. Only a few have really passed the test of time. Here is a summary of a modern proteomics analysis pipeline as implemented in the PAW pipeline:

1. Extract MS2 scan information from the RAW files (MSConvert [2, 3])
   * Convert to a format the search engine can read
1. Select an appropriate FASTA protein database
   * Avoid excessive peptide redundancy
   * Add decoys and/or contaminants if necessary
1. Correctly configure the search engine (Comet [4])
   * Choose appropriate parent and fragment ion tolerances
   * Avoid excessive post-translational modifications
1. Transform search scores into more sensitive functions
   * Combinations of scores can work better
   * Use accurate mass
1. Filter out the confident PSMs
   * Use target/decoy methods to control PSM false discovery rate (FDR)
1. Use basic or extended parsimony logic to infer proteins from peptides
   * Avoid single peptide per protein identifications
   * Avoid protein ranking functions
     * Protein FDR is a consequence of peptide ID accuracy
     * Protein error control is different than PSM error control
1. Add protein/gene annotations (optional)

Identification of the proteins present in a sample is almost never the end goal of a modern proteomics experiment. Estimating the relative expression levels of the proteins is often required. The above discussion and list has not mentioned anything about quantification. Quantitative processing is really more of a parallel set of analysis steps. Quantitative information can take many forms. There are label free approaches and stable isotope labeling approaches. There is no need to survey quantitative proteomics as the PAW pipeline only does tandem mass tag (TMT) labeling. The support for TMT is even more limited to high resolution instruments. In a similar inference process, protein expression values are inferred from quantitative data acquired in individual instrument scans.

Here is the way TMT quantification is supported:

1. Convert RAW files into compressed text format (MSConvert [2, 3])
1. Extract maximum intensity in narrow windows centered at reporter ion masses
   * Create mapping between MS2 and MS3 scan numbers
   * Save extracted data in text files
1. Read PAW results files to determine usable quantitative data
   * Protein inference determines unique peptide status
   * Peptide and PSM reports yield lists of scans for unique peptides
1. Total protein reporter ion intensities are computed from the scan lists
   * Minimum intensity level testing
   * Optional replacement of missing reporter ion intensities
1. Reporter ion intensities are added to the PAW results files
1. Pipeline supports multiple TMT labeling experiments
   * Internal reference scaling (IRS) normalization [18] can be performed

These lists demonstrate that even a basic processing pipeline will involve several steps. A robust pipeline will keep these step separate to allow greater flexibility and to allow inspection of the data in between steps for quality control.

## Pipeline Steps
The PAW pipeline is an actual pipeline. The steps are modular and separate by design.

![PAW workflow diagram(images/PAW_workflow.png)

Each step is listed below by the name of the Python script. Each step has a link to dedicated documentation pages.

* [`msconvert_GUI.py`](docs/msconvert_GUI.md) - Conversion of RAW files
  - also `MSConvert_GUI_guide.pptx`
* [`comet_GUI.py`](docs/comet_GUI.md) - Setting Comet parameters and running Comet
  - also `comet_GUI_guide.pptx`
* [`histo_GUI.py`](docs/histo_GUI.md) - Viewing score distributions and setting thresholds
  - also `histo_GUI_guide.pptx`
* [`PAW_results.py`](docs/PAW_results.md) - Protein inference and results report generation
  - `PAW_protein_grouper.py` - Extended parsimony protein grouping is run by above script
  - also `PAW_results_guide.pptx`
* [`add_TMT_intensities.py`](docs/add_TMT_intensities.md) - (optional) Adding reporter ion intensities to reports
  - also `add_TMT_intensities.pptx`
* [`pandas_TMT_IRS_norm.py`](docs/pandas_TMT_IRS_norm.md) - (optional) IRS normalization between TMT experiments

## Dependancies
The programs to (1) read the RAW instrument binary files and (2) perform the database search have to be installed separately and they only run on Windows (7 or 10) systems.

#### MSConvert
Thermo RAW file conversions are done using ProteoWizard [2, 3] and its `MSConvert.exe` program. ProteoWizard installation instructions can be found [**_>> HERE <<_**](docs/MSConvert.md).

#### Comet
Comet installation instructions can be found [**_>> HERE <<_**](docs/Comet.md). Comet is not strictly a Windows only program (Linux and Mac are possible); however, installation can be more involved and is not covered here.

#### Python 3
The Python scripts require Python 3 and some packages that are not part of the standard distribution. Use of a scientific Python distribution is recommended and installation and use of Anaconda is described [**_>> HERE <<_**](docs/Anaconda.md).

After RAW files have been processed and searches completed, the Python scripts can be run on Windows or Mac. A little care has to be paid to hard disc formats, however. RAW files and Comet results files will need to be on an NTFS formatted drive (the Windows file system). NTFS drives can be read by Macs, but they cannot be written to. You will either have to copy files to Mac formatted volumes or use the [Paragon NTFS for Mac](https://www.paragon-software.com/ufsdhome/store/ntfs-mac/) utility software.

---

## References
[1] Wilmarth, P.A., Riviere, M.A. and David, L.L., 2009. Techniques for accurate protein identification in shotgun proteomic studies of human, mouse, bovine, and chicken lenses. Journal of ocular biology, diseases, and informatics, 2(4), pp.223-234.

[2] Kessner, D., Chambers, M., Burke, R., Agus, D. and Mallick, P., 2008. ProteoWizard: open source software for rapid proteomics tools development. Bioinformatics, 24(21), pp.2534-2536.

[3] Chambers, M.C., Maclean, B., Burke, R., Amodei, D., Ruderman, D.L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J. and Hoff, K., 2012. A cross-platform toolkit for mass spectrometry and proteomics. Nature biotechnology, 30(10), p.918.

[4] Eng, J.K., Jahan, T.A. and Hoopmann, M.R., 2013. Comet: an open‐source MS/MS sequence database search tool. Proteomics, 13(1), pp.22-24

[5] https://github.com/pwilmart/fasta_utilities.git

[6] Elias, J.E. and Gygi, S.P., 2007. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. Nature methods, 4(3), p.207.

[7] Elias, J.E. and Gygi, S.P., 2010. Target-decoy search strategy for mass spectrometry-based proteomics. In Proteome bioinformatics (pp. 55-71). Humana Press.

[8] Hsieh, E.J., Hoopmann, M.R., MacLean, B. and MacCoss, M.J., 2009. Comparison of database search strategies for high precursor mass accuracy MS/MS data. Journal of proteome research, 9(2), pp.1138-1143.

[9] Keller, A., Nesvizhskii, A.I., Kolker, E. and Aebersold, R., 2002. Empirical statistical model to estimate the accuracy of peptide identifications made by MS/MS and database search. Analytical chemistry, 74(20), pp.5383-5392.

[10] Searle, B.C., 2010. Scaffold: a bioinformatic tool for validating MS/MS‐based proteomic studies. Proteomics, 10(6), pp.1265-1269.

[11] Nesvizhskii, A.I. and Aebersold, R., 2005. Interpretation of shotgun proteomic data the protein inference problem. Molecular & cellular proteomics, 4(10), pp.1419-1440.

[12] Madhira, R., 2016. The Effects of Parsimony Logic and Extended Parsimony Clustering on Protein Identification and Quantification in Shotgun Proteomics. https://digitalcollections.ohsu.edu/concern/etds/c534fp149

[13] Deutsch, E.W., Sun, Z., Campbell, D.S., Binz, P.A., Farrah, T., Shteynberg, D., Mendoza, L., Omenn, G.S. and Moritz, R.L., 2016. Tiered human integrated sequence search databases for shotgun proteomics. Journal of proteome research, 15(11), pp.4091-4100.

[14] Eng, J.K., McCormack, A.L. and Yates, J.R., 1994. An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database. Journal of the American Society for Mass Spectrometry, 5(11), pp.976-989.

[15] Tabb, D.L., McDonald, W.H. and Yates, J.R., 2002. DTASelect and Contrast: tools for assembling and comparing protein identifications from shotgun proteomics. Journal of proteome research, 1(1), pp.21-26.

[16] Käll, L., Canterbury, J.D., Weston, J., Noble, W.S. and MacCoss, M.J., 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nature methods, 4(11), p.923.

[17] Carr, S., Aebersold, R., Baldwin, M., Burlingame, A.L., Clauser, K. and Nesvizhskii, A., 2004. The need for guidelines in publication of peptide and protein identification data working group on publication guidelines for peptide and protein identification data.

[18] Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

[to top](#paw_pipeline)
