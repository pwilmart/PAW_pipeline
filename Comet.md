# Comet installation and setup

Comet is an open source search engine written by Jimmy Eng at the University of Washington. The project is [hosted at SourceForge](http://comet-ms.sourceforge.net/) and the downloads page is located [**_>> HERE <<_**](https://sourceforge.net/projects/comet-ms/files/). Jimmy authored one of the first and arguably best search engines, SEQUEST [2], and is a principle scientist at the [UW Proteomics Resource](http://proteomicsresource.washington.edu/).

Comet binaries download as `.ZIP` archives and will need to be located in an appropriate folder and the `.ZIP` archive uncompressed:


![Comet_download](images/Comet_download.png)

![Comet_folder](images/Comet_folder.png)

Comet undergoes occasional updates and bug fixes. One way to isolate the PAW Python scripts from the desired Comet version, is to have a BATCH file point to the desired Comet executable file. If a `.BAT` file is in a folder located in the system search path, then the base name of the `.BAT` file can be called by the script (and stay the same), and execution of the current version of the software can be maintained by editing the program path in the `.BAT` file.

![Comet_batch](images/Comet_batch.png)

<br>

Here are the contents of the `Comet.bat` file:
```
C:\Users\pwilmart\Proteomics_tools\Comet\comet_2017014\comet.2017014.win64.exe %*
```

<br>

The path in the `Comet.bat` file should match the latest installed version of `Comet`. The `%*` at the end of the line passes any batch command line arguments on to the `comet.2017014.win64.exe` program. The `commands_scripts` folder contains batch files and some other utilities and has been added to the system search path. Any batch files located here can be executed from any Windows Command or PowerShell window, or be called from scripting languages like Python. For instructions on how to add a folder to the system path in Windows 10 see [**THIS LINK**](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho).

When everything is in place, you can test the the `Comet.bat` file can be found by Windows. Open a  Command or PowerShell window (SHIFT-RIGHT-CLICK on a folder to open a shell window with the default location set to the folder) and type `Comet` at the prompt. You will get a usage message if everything is working correctly.

![Comet_command](images/Comet_command.png)

<br>

## References

[1] Eng, J.K., Jahan, T.A. and Hoopmann, M.R., 2013. Comet: an open‐source MS/MS sequence database search tool. Proteomics, 13(1), pp.22-24

[2] Eng, J.K., McCormack, A.L. and Yates, J.R., 1994. An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database. Journal of the American Society for Mass Spectrometry, 5(11), pp.976-989.