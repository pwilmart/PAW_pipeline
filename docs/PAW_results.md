# PAW_results.py script

This script takes the filtered top hit files and makes a parsimonious list of identified proteins following the basic logic outlined here:

> Nesvizhskii, A.I. and Aebersold, R., 2005. Interpretation of shotgun proteomic data the protein inference problem. Molecular & cellular proteomics, 4(10), pp.1419-1440.

Multiple fractions per biological sample and multiple biological samples per experiment are supported. Protein inference is performed experiment-wide to provide robust results. Optional [extended parsimony logic](https://digitalcollections.ohsu.edu/record/3031?ln=en&p=Madhira&v=pdf) is also implemented to group highly homologous proteins and maximize the quantitative information content of peptides. Final protein and peptide reports are complete, carefully curated tab-delimited text files.

### Step-by-step Example, continued

This is analysis of a public dataset (PRIDE [PXD002875](https://www.ebi.ac.uk/pride/archive/projects/PXD002875)) from Paulo, O'Connell, Gaun, and Gygi:

> Paulo, J.A., O’Connell, J.D., Gaun, A. and Gygi, S.P., 2015. Proteome-wide quantitative multiplexed profiling of protein expression: carbon-source dependency in Saccharomyces cerevisiae. Molecular biology of the cell, 26(22), pp.4063-4074.

There were 24 RAW files of yeast grown in three different carbon sources. It was a 3x3 (9-plex) TMT experiment done with the SPS MS3 (MultiNotch) method.

**Processing flow:**
- read filtered top-hit summary results_files
- define file names to samples mapping
- perform basic parsimonious protein inference
- write results files
- perform extended parsimony protein grouping
- write grouped protein results files

---

![start anaconda prompt](../images/PAW_results/01_anaconda.png)

Start an Anaconda Prompt window

---

![start IDLE](../images/PAW_results/02_idle.png)

Type the command `idle` and press return.

---

![open script](../images/PAW_results/03_open-script.png)

This will launch the IDLE programming environment. From the IDLE console window, select the `Open...` command from the `File` menu.

---

![select script](../images/PAW_results/04_select-script.png)

Select the `PAW_results.py` script from the location where you have the PAW pipeline scripts from the file browser dialog box. Then click the `Open` button in the lower right corner.

---

![parameters](../images/PAW_results/05_parameters.png)

There are a couple of parameters that you might want to occasionally change. Two peptides per protein is a **HIGHLY recommended** setting. Sometimes you may need single peptide per protein IDs. The PAW pipeline does not have any explicit **protein** FDR algorithms. There is a one-to-one relationship between incorrect PSMs and incorrect proteins when using a one peptide per protein setting. Sometimes semi-tryptic or non-tryptic peptides are important. The `minimum_ntt_per_peptide` parameter would need adjustment in those cases.

---

![run script](../images/PAW_results/06_run-script.png)

Run the script (`Run Module` command from the `Run` menu).

---

![select folder](../images/PAW_results/07_select-folder.png)

The IDLE console windows indicates that a folder dialog box has been popped up for `filtered_files` folder selection.

---

![folder dialog](../images/PAW_results/08-folder-dialog.png)

Browse to the `filtered_files` folder inside your project folder. Select the folder and then click the `Select Folder` button.

---

![CLI help](../images/PAW_results/09_CLI-help.png)

All of the filtered top hit summary files are parsed and a simple command line interpreter (CLI) is started. This allows assigning the RAW files to the biological samples. Many experiment have multiple biological samples. Shotgun analysis of biological samples often requires multiple dimensions of separation, so several RAW files might correspond to single biological samples. We have the “Help” information in red above. The pattern matching is like most file system [glob pattern](https://en.wikipedia.org/wiki/Glob_(programming)) matching.

---

![CLI list](../images/PAW_results/10_CLI-list.png)

A good command to do first is to `List` the files. We have the 24 fractions for this sample. Since this is a 9-plex TMT experiment, the 9 biological samples are hidden. The 24 fractions are really just one effective biological sample.

---

![CLI pattern](../images/PAW_results/11_CLI-pattern.png)

We will use the `Pattern` command to select all of the files so we can give them a nice name. The pattern `m*` will select all files starting with an “m” character.

---

![pattern matches](../images/PAW_results/12_pattern-matches.png)

The 24 matching files are echoed and a sample name is asked for. We will use “CarbonSources” for the effective sample name.

---

![CLI show](../images/PAW_results/13_CLI-show.png)

If we want to double check what has been assigned to samples so far, we can type the `Show` command. It tells us that we have one sample with 24 files. Collecting files to assign to samples is usually an iterative process.

---

![CLI done](../images/PAW_results/14_CLI-done.png)

We  can make sure that we have no leftover files (files not yet assigned to samples) by typing the `List` command again. It shows that we do not have any remaining files. We type the `Done` command to exit the command processor and get back to the protein inference algorithm.

---

![load results](../images/PAW_results/15_load-results.png)

The files to samples mapping(s) have been defined, so we can read in the filtered results files. After the files have been read in, we get a PSM FDR report. It is broken down by charge state and by modification state.

---

![redundant proteins](../images/PAW_results/16_redundants.png)

After the FDR report, the protein inference starts by making peptide sets for each protein that had any PSMs that passed the thresholds. Any proteins having peptide sets with fewer distinct peptides than the minimum cutoff of 2 are removed. The next step is to group proteins having indistinguishable peptide sets (also called redundant proteins).

---

![subset proteins](../images/PAW_results/17_subsets.png)

Grouping redundant proteins decreased the number of proteins from 5150 to 5028. Next, we do the main parsimony step and remove proteins whose peptide sets are subsets of proteins having larger peptide sets.

---

![half done](../images/PAW_results/18_half-done.png)

After subset removal, we have 4945 protein groups. Yeast has a simple genome, so protein grouping does not alter the protein numbers as much as we see for species like human and mouse. The protein grouping algorithm runs automatically.

---

![protein grouping](../images/PAW_results/19_protein-grouping.png)

The protein grouping algorithm is an extension of the basic parsimony concepts. If two proteins have nearly identical peptides sets, they will be grouped. If one protein has a peptide set that is almost a subset of a protein with a larger peptide set, the subset protein will be absorbed. The test conditions are weight tests. The weights are the total spectral counts. When the **weight** of one protein’s unique and shared peptides are much larger than another protein’s unique **weight**, protein grouping will happen. The grouping is iterative and continues until the final number of groups is no longer changing. Grouped proteins get their peptide’s shared and unique status redefined after grouping. This can have some effects on quantitative results for highly homologous protein families.

---

![protein families](../images/PAW_results/20_protein-families.png)

The protein grouping console output is usually really interesting to look at. All details during the iterative steps are logged. After the final set of protein families has been determined, the family memberships are logged. The protein descriptions for the grouped proteins indicated similar proteins. Housekeeping genes are frequently grouped.

---

![final summary](../images/PAW_results/21_final-summary.png)

Summary statistics are logged at the end.

More details on the concepts behind the extended parsimony grouping can be found at this [Cascadia 2013 talk](https://github.com/pwilmart/Cascadia_2013) repository.

You can watch an example of the algorithm in action in [this video link](https://youtu.be/IAE59QrQYos). It is just over 2 minutes long. The data is from [this repository](https://github.com/pwilmart/MaxQuant_and_PAW/tree/master/PAW_results) and the `PAW_protein_grouper.log` file has all of the console output. The samples are from cultured mouse bone marrow cells.

---

![results folder](../images/PAW_results/22_results-folder.png)

We have a new `results_files` folder added to our project folder.

---

![results files](../images/PAW_results/23_results-files.png)

That folder has our tab-delimited text summary files. There will be a detailed peptide file for each biological sample (we have only one here: `CarbonSource_peptide_resuls_8.txt`). We have the redundant protein and peptide reports generated before protein grouping. We have the non-redundant reports produced by the protein grouping step, and we have log files that duplicate most of the console output.

> There is a `CarbonSources_results_files` folder inside the `docs` folder in the repository that has the results files created during this example. Your results may not be 100% identical because your thresholds set during the histo_GUI step may not have been exactly the same as mine. Results should be very similar, however.

Optional next step: [Add the TMT reporter ions to the reports](add_TMT_intensities.md)
