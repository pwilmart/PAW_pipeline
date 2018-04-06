# PAW_pipeline
A Comet-based, best practices proteomics pipeline.

* [Introduction](#introduction)
* [Background](#background)
* [Best practices](#best-practices)
* [Dependancies](#dependancies)
* [References](#references)

## Introduction
The PAW pipeline [1] consists of open source Python scripts to process Thermo mass spectrometer RAW files in a best-practices shotgun proteomics analysis pipeline. The open source ProteoWizard toolkit [2, 3] is used to convert scans in the RAW files into more easily readable formats. High-resolution isobaric labeling data (10- or 11-plex TMT reagents) can also be handled. Scripts simplify configuring and running MSConvert, convert the MS2 scans into appropriate file formats for database search, and extract the TMT reporter ion peak heights.  

The open source Comet database search engine [4] is used to identify peptides. Scripts simplify setting relevant search parameters, launching the Comet searches, and creating top-hit summary files. Other scripts [5] assist in downloading high-quality, up-to-date protein FASTA databases and adding decoy sequences and common contaminants.

Comet results are accurately filtered to desired false discovery rates (FDRs) using the target decoy method [6, 7] and proper use of accurate mass [8]. A GUI script generates accurate-mass-conditioned score histograms and FDR tables. The data are subclassed by peptide charge state and modification state. The scores are a linear discriminant function transformation similar to PeptideProphet [9] and Scaffold [10]. Score thresholds are set interactively in each data subclass. Inspection of score histograms provides integrated quality control.

Protein inference follows basic parsimony rules [11]. Multiple fractions per biological sample and multiple biological samples per experiment are supported. Protein inference is performed experiment-wide to provide robust results. Optional extended parsimony logic [12] is also implemented to group highly homologous proteins and maximize the quantitative information content of peptides. Reporter ion quantitative information in TMT experiments are aggregated to protein expression totals. Final protein and peptide reports are complete, carefully curated tab-delimited text files.

## Background
Analyzing bottom-up (a.k.a. shotgun) proteomics data can be a complicated, lengthy ordeal. The processing pipelines have many steps and each step can have several choices. There are few practical guidelines for making better choices to navigate this minefield. This results in far too great a diversity of pipelines in use, and most will have less-than-ideal performance.

This state of affairs is illustrated by one of the first steps of selecting a FASTA protein database file for use by the search engine. This is a more complicated topic than might be guessed [12]. There are often several protein database choices for many commonly studied organisms (and even some genomes that have not yet been sequenced). The protein database complexity and apparent size can vary considerably among the choices, although larger is probably not better [13]. Even integrated database sources like [UniProt](http://www.uniprot.org/) offer several versions of protein databases (Swiss-Prot, TrEMBL, canonical, canonical with isoforms, etc.), and selecting among them requires some expertise and use of [dedicated tools](https://github.com/Delan-Huang/Reference_Proteome_Manager.git) [5].  

Picking a search engine program to identify likely peptide sequences associated with tandem mass spectra, known as peptide spectrum matches (PSMs), can be even more challenging. There are commercial products like Mascot, Proteome Discoverer, and Byonic that can be quite expensive. Alternatively, there are widely used freely available options such as X!Tandem and MaxQuant. There are also many less commonly used open source options, and opinions run high on the relative merits of the various search engines. Comet [4], written and maintained by Jimmy Eng, is an open source version of SEQUEST [14]. Comet is free, fast, sensitive, and a little simpler to configure than other options; this is what the PAW pipeline [1] uses. PAW are my initials, but you can make up something like Proteomic Analysis Workflow if you prefer.

The basic method of identifying peptides from tandem mass spectra (MS2 scans) has not really changed too much in 25 years [14]. It basically involves noise filtering and normalizing an MS2 scan, then comparing that to mass-filtered candidate peptide spectra from a theoretical digest of a protein database. Despite this clever leveraging of genomic sequencing (the protein sequences), the challenge has always been in deciding if a particular PSM is correct or not. Early heuristic approaches [15] led to basic classifiers [9], and even machine learning methods [16]. A big step forward came when decoy databases were used to eavesdrop on random noise scores, as popularized by the Gygi lab [6, 7]. More recent instrumental advances have added high-resolution and accurate mass to the equation [8].

With more robust statistical methods for controlling PSM errors, confident lists of identified **peptides** present in digests of complex protein mixtures can be reliably determined. The next issue is to determine a confident list of the **proteins** present in the samples from the identified peptides. This is known as the protein inference problem [11] and it persists despite significant advances in genomic sequencing (maybe that has made the problem even worse). Guidelines for parsimonious protein identifications [17] have been widely used for many years in response to many early inflated proteome claims. The basic parsimony rules are outlined in [11], but may need to be extended [12] for the large datasets now routinely being generated.

## Best Practices
As proteomics has matured, there have been many analysis ideas that have come and gone. Only a few have really passed the test of time. Here is a summary of a modern proteomics analysis pipeline as implemented in the PAW pipeline:

1. Extract MS2 scan information from the RAW files (MSConvert [2, 3])
   * Convert to a format the search engine can read
1. Select an appropriate FASTA protein database
   * Avoid excessive peptide redundancy
   * Add decoys and/or contaminants if necessary
1. Correctly configure the search engine (Comet [4])
   * Choose parent and fragment ion tolerances wisely
   * Avoid excessive post-translational modifications
1. Transform search scores into more sensitive functions
   * Combinations of scores can work better
   * Use accurate mass wisely
1. Filter out the confident PSMs
   * Use target/decoy methods to control PSM false discovery rate (FDR)
1. Use basic or extended parsimony logic to infer proteins from peptides
   * Avoid single peptide per protein identifications
   * Avoid using protein ranking functions
     * Protein FDR is a consequence of peptide ID accuracy
     * Protein error control is different than PSM error control
1. Add protein/gene annotations (optional)

Identification of the proteins present in a sample is almost never the goal of a modern proteomics experiment. Estimating the relative expression levels of the proteins is often required. The above discussion and list has not mentioned anything about quantification. Quantitative processing is really more of a parallel set of analysis steps. Quantitative information can take many forms. There are label free approaches and stable isotope labeling approaches. There is no need to survey quantitative proteomics as the PAW pipeline only does tandem mass tag (TMT) labeling. The support for TMT is even more limited to high resolution instruments and synchronous precursor selection (SPS) MS3 scan  reporter ions. In a similar inference process, protein expression values are inferred from quantitative data acquired in individual instrument scans.

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

These lists demonstrate that even a basic processing pipeline will involve several steps. A robust pipeline will keep these step separate to allow greater flexibility and to allow inspection of the data in between steps to quality control. The PAW software was designed to **not** try and cram 10 pounds of pipeline into a 5 pound _black box_.

## Dependancies
The programs to read the RAW instrument binary files and perform the database search have to be installed separately and only run on Windows systems.

Thermo RAW file conversions are done using ProteoWizard [2, 3] and its `MSConvert.exe` program. ProteoWizard installation instructions can be found [**_>> HERE <<_**](MSConvert.md).

Comet installation instructions can be found [**_>> HERE <<_**](Comet.md). Comet is not strictly a Windows only program (Linux and Mac are possible); however, installation can be more involved and is not covered here.

The Python scripts require Python 3 and some packages that are not part of the standard distribution. Use of a scientific Python distribution is recommended and installation and use of Anaconda is described [**_>> HERE <<_**](Anaconda.md). After RAW files have been processed and searches completed, the Python scripts can be run on Windows or Mac. A little care has to be paid to hard disc formats, however. RAW files and Comet results files will need to be on an NTFS formatted drive (the Windows file system). NTFS drives can be read by Macs, but they cannot be written to. You will either have to copy files to Mac formatted volumes or use the [Paragon NTFS for Mac](https://www.paragon-software.com/ufsdhome/store/ntfs-mac/) utility software.

<br>

## References
[1] Wilmarth, P.A., Riviere, M.A. and David, L.L., 2009. Techniques for accurate protein identification in shotgun proteomic studies of human, mouse, bovine, and chicken lenses. Journal of ocular biology, diseases, and informatics, 2(4), pp.223-234.

[2] Kessner, D., Chambers, M., Burke, R., Agus, D. and Mallick, P., 2008. ProteoWizard: open source software for rapid proteomics tools development. Bioinformatics, 24(21), pp.2534-2536.

[3] Chambers, M.C., Maclean, B., Burke, R., Amodei, D., Ruderman, D.L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J. and Hoff, K., 2012. A cross-platform toolkit for mass spectrometry and proteomics. Nature biotechnology, 30(10), p.918.

[4] Eng, J.K., Jahan, T.A. and Hoopmann, M.R., 2013. Comet: an open‐source MS/MS sequence database search tool. Proteomics, 13(1), pp.22-24

[5] https://github.com/Delan-Huang/Reference_Proteome_Manager.git

[6] Elias, J.E. and Gygi, S.P., 2007. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. Nature methods, 4(3), p.207.

[7] Elias, J.E. and Gygi, S.P., 2010. Target-decoy search strategy for mass spectrometry-based proteomics. In Proteome bioinformatics (pp. 55-71). Humana Press.

[8] Hsieh, E.J., Hoopmann, M.R., MacLean, B. and MacCoss, M.J., 2009. Comparison of database search strategies for high precursor mass accuracy MS/MS data. Journal of proteome research, 9(2), pp.1138-1143.

[9] Keller, A., Nesvizhskii, A.I., Kolker, E. and Aebersold, R., 2002. Empirical statistical model to estimate the accuracy of peptide identifications made by MS/MS and database search. Analytical chemistry, 74(20), pp.5383-5392.

[10] Searle, B.C., 2010. Scaffold: a bioinformatic tool for validating MS/MS‐based proteomic studies. Proteomics, 10(6), pp.1265-1269.

[11] Nesvizhskii, A.I. and Aebersold, R., 2005. Interpretation of shotgun proteomic data the protein inference problem. Molecular & cellular proteomics, 4(10), pp.1419-1440.

[12] Madhira, R., 2016. The Effects of Parsimony Logic and Extended Parsimony Clustering on Protein Identification and Quantification in Shotgun Proteomics. https://digitalcommons.ohsu.edu/etd/3855/

[13] Deutsch, E.W., Sun, Z., Campbell, D.S., Binz, P.A., Farrah, T., Shteynberg, D., Mendoza, L., Omenn, G.S. and Moritz, R.L., 2016. Tiered human integrated sequence search databases for shotgun proteomics. Journal of proteome research, 15(11), pp.4091-4100.

[14] Eng, J.K., McCormack, A.L. and Yates, J.R., 1994. An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database. Journal of the American Society for Mass Spectrometry, 5(11), pp.976-989.

[15] Tabb, D.L., McDonald, W.H. and Yates, J.R., 2002. DTASelect and Contrast: tools for assembling and comparing protein identifications from shotgun proteomics. Journal of proteome research, 1(1), pp.21-26.

[16] Käll, L., Canterbury, J.D., Weston, J., Noble, W.S. and MacCoss, M.J., 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nature methods, 4(11), p.923.

[17] Carr, S., Aebersold, R., Baldwin, M., Burlingame, A.L., Clauser, K. and Nesvizhskii, A., 2004. The need for guidelines in publication of peptide and protein identification data working group on publication guidelines for peptide and protein identification data.

[18] Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

[to top](#paw_pipeline)
