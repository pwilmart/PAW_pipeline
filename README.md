# PAW_pipeline
A Comet-based, best practices proteomics pipeline.

* [Background](#background)
* [References](#references)



## Background
Bottom-up (a.k.a. shotgun) proteomics data analysis can be a complicated, lengthy ordeal. The processing pipelines have many steps and each step can have several choices. There are few practical guidelines for making better choices to navigate this minefield. This results in far too great a diversity of pipelines in use, and many may have less than ideal performance.

This state of affairs is perfectly illustrated by the first step of selecting a FASTA protein database file for use by the search engine. This is a more complicated topic than might be guessed [1]. There are several choices for many organisms and even some genomes that have not yet been sequenced. The protein database complexity and apparent size can vary considerably, although larger is probably not better [2]. Even single database sources like [UniProt](www.uniprot.org) offer several versions of protein databases and selecting among them can require [dedicated tools](https://github.com/Delan-Huang/Reference_Proteome_Manager.git) [3].  

Picking a search engine program to identify likely peptide sequences associated with tandem mass spectra, known as peptide spectrum matches (PSMs), can be even more challenging. There are commercial products like Mascot, Proteome Discoverer, and Byonic that can be quite expensive. Alternatively, there are widely used freely available options such as X!Tandem and MaxQuant. There are also many less commonly used open source options and opinions run high on the relative merits of the various search engines. Comet [4], written and maintained by Jimmy Eng, is an open source version of SEQUEST. Comet is free, fast, sensitive, and a little simpler to configure that other options; this is what the PAW pipeline [5] uses (PAW are my initials but you can make up something like Proteomic Analysis Workflow if you prefer).

The basics of identifying peptides from tandem mass spectra (MS2 scans) has not really changed conceptually too much in 25 years [6]. It basically involves noise filtering and normalizing an MS2 scan, then comparing that to mass-filtered candidate peptide spectra from a theoretical digest of a protein database. Despite this clever leveraging of genomic sequencing (the protein sequences), the challenge has always been in deciding if a particular PSM is correct or not. Early heuristic approaches [7] led to basic classifiers [8], and even machine learning methods [9]. A key step forward was use of decoy databases to eavesdrop on random noise scores that was popularized by the Gygi group [10, 11]. More recent instrumental advances have added high-resolution and accurate mass to the equation [12].

With more robust statistical methods for controlling PSM errors, confident lists of identified **peptides** present in digests of complex protein mixtures can be reliably determined. The next issue is to determine a confident list of the **proteins** present in the samples from the identified peptides. This is known as the protein inference problem [13] and it persists despite significant advances in genomic sequencing (one could argue that the problem is even worse). Guidelines for parsimonious protein identifications [14] have been widely used for many years in response to many early inflated proteome claims. The basic parsimony rules are outlined in [13], but may need to be extended [1] in the context of the large datasets now routinely being generated.



<br>

## References:
[1] https://digitalcommons.ohsu.edu/etd/3855/

[2] Deutsch, E.W., Sun, Z., Campbell, D.S., Binz, P.A., Farrah, T., Shteynberg, D., Mendoza, L., Omenn, G.S. and Moritz, R.L., 2016. Tiered human integrated sequence search databases for shotgun proteomics. Journal of proteome research, 15(11), pp.4091-4100.

[3] https://github.com/Delan-Huang/Reference_Proteome_Manager.git

[4] Eng, J.K., Jahan, T.A. and Hoopmann, M.R., 2013. Comet: an open‐source MS/MS sequence database search tool. Proteomics, 13(1), pp.22-24

[5] Wilmarth, P.A., Riviere, M.A. and David, L.L., 2009. Techniques for accurate protein identification in shotgun proteomic studies of human, mouse, bovine, and chicken lenses. Journal of ocular biology, diseases, and informatics, 2(4), pp.223-234.

[6] Eng, J.K., McCormack, A.L. and Yates, J.R., 1994. An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database. Journal of the American Society for Mass Spectrometry, 5(11), pp.976-989.

[7] Tabb, D.L., McDonald, W.H. and Yates, J.R., 2002. DTASelect and Contrast: tools for assembling and comparing protein identifications from shotgun proteomics. Journal of proteome research, 1(1), pp.21-26.

[8] Keller, A., Nesvizhskii, A.I., Kolker, E. and Aebersold, R., 2002. Empirical statistical model to estimate the accuracy of peptide identifications made by MS/MS and database search. Analytical chemistry, 74(20), pp.5383-5392.

[9] Käll, L., Canterbury, J.D., Weston, J., Noble, W.S. and MacCoss, M.J., 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nature methods, 4(11), p.923.

[10] Elias, J.E. and Gygi, S.P., 2007. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. Nature methods, 4(3), p.207.

[11] Elias, J.E. and Gygi, S.P., 2010. Target-decoy search strategy for mass spectrometry-based proteomics. In Proteome bioinformatics (pp. 55-71). Humana Press.

[12] Hsieh, E.J., Hoopmann, M.R., MacLean, B. and MacCoss, M.J., 2009. Comparison of database search strategies for high precursor mass accuracy MS/MS data. Journal of proteome research, 9(2), pp.1138-1143.

[13] Nesvizhskii, A.I. and Aebersold, R., 2005. Interpretation of shotgun proteomic data the protein inference problem. Molecular & cellular proteomics, 4(10), pp.1419-1440.

[14] Carr, S., Aebersold, R., Baldwin, M., Burlingame, A.L., Clauser, K. and Nesvizhskii, A., 2004. The need for guidelines in publication of peptide and protein identification data working group on publication guidelines for peptide and protein identification data.

Kessner, D., Chambers, M., Burke, R., Agus, D. and Mallick, P., 2008. ProteoWizard: open source software for rapid proteomics tools development. Bioinformatics, 24(21), pp.2534-2536.

Chambers, M.C., Maclean, B., Burke, R., Amodei, D., Ruderman, D.L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J. and Hoff, K., 2012. A cross-platform toolkit for mass spectrometry and proteomics. Nature biotechnology, 30(10), p.918.

[to top](#pawpipeline)
