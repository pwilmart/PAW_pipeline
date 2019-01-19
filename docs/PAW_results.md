

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

There are a couple of parameters that you might want to occasionally change. Two peptides per protein is a **HIGHLY recommended** setting. Sometimes you may need single peptide per protein IDs. The PAW pipeline does not have any explicit protein FDR algorithms. There is a one-to-one relationship between incorrect PSMs and incorrect proteins when using a one peptide per protein setting. Sometimes semi-tryptic or non-tryptic peptides are important. The “minimum_ntt_per_peptide” parameter would need adjustment in those cases.

---

![run script](../images/PAW_results/06_run-script.png)

Run the script (`Run Module` command from the `Run` menu).

---

![select folder](../images/PAW_results/07_select-folder.png)

The IDLE console windows indicates that a folder dialog box has been popped up for “filtered_files” folder selection.

---

![folder dialog](../images/PAW_results/08-folder-dialog.png)

Browse to the `filtered_files` folder inside your project folder. Select the folder and then click the `Select Folder` button.

---

![CLI help](../images/PAW_results/09_CLI-help.png)

All of the filtered top hit summary files are parsed and a simple command interpreter is started. This allows assigning the RAW files to the biological samples. Many experiment have multiple biological samples. Shotgun analysis of biological samples often requires multiple dimensions of separation, so several RAW files might correspond to single biological samples. We have the “Help” information on this slide. The pattern matching is like most file system [glob pattern](https://en.wikipedia.org/wiki/Glob_(programming)) matching.


---

![CLI list](../images/PAW_results/10_CLI-list.png)

---

![CLI pattern](../images/PAW_results/11_CLI-pattern.png)

---

![pattern matches](../images/PAW_results/12_pattern-matches.png)

---

![CLI show](../images/PAW_results/13_CLI-show.png)

---

![CLI done](../images/PAW_results/14_CLI-done.png)

---

![load results](../images/PAW_results/15_load-results.png)

---

![redundant proteins](../images/PAW_results/16_redundants.png)

---

![subset proteins](../images/PAW_results/17_subsets.png)

---

![half done](../images/PAW_results/18_half-done.png)

---

![protein grouping](../images/PAW_results/19_protein-grouping.png)

---

![protein families](../images/PAW_results/20_protein-families.png)

---

![final summary](../images/PAW_results/21_final-summary.png)

---

![results folder](../images/PAW_results/22_results-folder.png)

---

![results files](../images/PAW_results/23_results-files.png)

---


```
(Auto Each Help List Pattern Reset Show Done) ? help

Commands: Auto, Each, Help, List, Pattern, Reset, Show, Quit.
(commands are not case sensitive and first letter is sufficient.)

...Auto: automatically parse the sample names from the filenames,
...Each: treat each filename as a separate sample,
...Help: prints this message,
...List: lists all of the (remaining) filename(s),
...Pattern pattern: glob-style search pattern to get a subset of filenames
      ("*" or "*.*" is all files,
       "*.txt" is all files that have a "txt" extension,
       "a*" is all files that start with the letter "a",
       "*string*" is all files that contain "string", etc.),
...Reset: start over,
...Show: print the current samples and associated files [verbose],
...Done: exit from the command loop and return to processing.
```
