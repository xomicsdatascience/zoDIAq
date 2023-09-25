# csdev_tdd_rewrite

This is the dev branch for the zodiaq 2.0 rewrite. In summary, we are enhancing the modularity and testability of zodiaq.

# zoDIAq
## Introduction

*Uses an MIT License* 

zoDIAq (Cosine Similarity Optimization for DIA qualitative and quantitative analysis, pronounced "Zodiac") is a software tool 
for the analysis of direct infusion shotgun proteome analysis (DISPA) data and data independent acquisition (DIA) data.

[Click here for the paper published in *Analytical Chemistry*.](https://doi.org/10.1021/acs.analchem.1c02021)

![Paper Title Header](https://github.com/xomicsdatascience/zoDIAq/blob/main/img/zodiaq-paper.png)

## Instructions
### Installation

All installation instructions are done from the command line. Download and use of zoDIAq requires access to the `git` and `pip` packages. For system-specific instructions on using the command line or installing these packages, see the [wiki page.](https://github.com/CCranney/zoDIAq/wiki)

* [Windows install instructions click here](https://github.com/CCranney/zoDIAq/wiki/Install-Instructions-for-Windows)
* [Mac install instructions click here](https://github.com/CCranney/zoDIAq/wiki/Install-Instructions-for-Mac)

### Accessing the GUI

From the command line enter `zodiaq gui` to start the GUI. Aside from the MGF library prep all instructions in this README file pertain to use of the GUI. 
Instructions for use from the command line can be accessed by entering `zodiaq -h`.

### Peptide/Protein Identification

![Peptide/Protein Identification GUI Picture](https://github.com/CCranney/zoDIAq/blob/master/Python%20Extras/ID_pic.png)

Files: 
1. DIA Data File: This is a required field. Choose at least one data file from a DISPA or DIA data run.
2. Library File: This is a required field. Choose a reference library file that is in TraML (.tsv or .csv) or MGF format. 
It should go without saying, but zoDIAq treats .tsv files as tab-deliminated and .csv as comma-deliminated. 
Some pan-human .csv files are tab-deliminated and therefore need to be adapted accordingly.
3. Outfile Directory: This is a required field. Choose a folder that the output files should go into.

Settings:

Note that it is recommended that the default value is used for each of these settings.

1. Initial Fragment Mass Tolerance (in PPM): Default value is 30. 
If you're not sure what would make an ideal setting, go ahead and leave this blank, but check the histogram box below. 
If the resulting histogram doesn't have a normal distribution (a peak), consider using a wider tolerance.
2. Correction: Default is checked. This value resets the PPM tolerance range based on an initial scan. It reduces the liklihood of identifying false positives.
3. Corrective Standard Deviation: Only available if the correction box is checked. Default is to customize the corrected tolerance to the distribution of the histogram, excluding noise around the peak. 
Entering a value here instead sets the tolerance to a standard deviation of the distribution.
4. Create Histogram: Only available if the correction box is checked. Default is checked. 
Generates a histogram to visualize the corrective PPM tolerance used in the corrected analysis.
5. Number of Target Peptides per Protein: Only available if the protein inference box is checked. 
Default is to use one target peptide to represent the protein in the output.
6. Maximum Number of Query Spectra to Pool: Default is to pool all query spectra. 
If you DISPA or DIA data file is particularly large and has a high number of scans with overlapping m/z windows, 
this setting restricts the number of scans that are pooled together for analysis to avoid exceeding your computer's memory capacity. 
Setting this value can slow down the program.
7. Permit Heavy Targets in Re-Analysis File: Default is checked. 
If checked, files generated for targeted re-analysis will include targets for heavy peptides (includes heavy lysine and arginine).

### SILAC Quantification

![SILAC Quantification GUI Picture](https://github.com/CCranney/zoDIAq/blob/master/Python%20Extras/quant_pic.png)


Files:
1. DIA Data File: This is a required field. Choose at least one data file from a DISPA or DIA targetted re-analysis data run.
2. Library File: This is a required field. Choose a reference library file that is in TraML (.tsv or .csv) or MGF format. 
It should go without saying, but zoDIAq treats .tsv files as tab-deliminated and .csv as comma-deliminated. 
Some pan-human .csv files are tab-deliminated and therefore need to be adapted accordingly.
3. Outfile Directory: This is a required field. Choose a folder that the output files should go into.
4. zoDIAq ID Output File: This is a required field. 
Choose the required output from the Peptide/Protein Identification portion. The file should end with "allCVs."

Settings:
1. Initial Fragment Mass Tolerance (in PPM): Default value is 30. 
If you're not sure what would make an ideal setting, go ahead and leave this blank, but check the histogram box below. 
If the resulting histogram doesn't have a normal distribution (a peak), consider using a wider tolerance.
2. Correction: Default is checked. This value resets the PPM tolerance range based on an initial scan. It reduces the liklihood of identifying false positives.
3. Corrective Standard Deviation: Only available if the correction box is checked. Default is to customize the corrected tolerance to the distribution of the histogram, excluding noise around the peak. 
Entering a value here instead sets the tolerance to a standard deviation of the distribution.
4. Create Histogram: Only available if the correction box is checked. Default is checked. 
Generates a histogram to visualize the corrective PPM tolerance used in the corrected analysis.
5. Number of Max Peaks per Library Spectra: Default is to use all library peaks. 
Changing this setting will sort the library peaks by intensity and only include the top __ values, as provided in this setting.
It is recommended that this value be set if the number of min peak matches setting is not set to the default.
6. Number of Min Peak Matches Required: Default is to require matching to one of the top 3 most intense peaks. 
This is the only setting that then includes all other matched peaks after the initial match - setting a value restricts the final values to the number set here.
7. Ratio Selection Method: Default is median. 
When determining the ratio between a library and query spectrum, ratios are calculated for each matching peak. 
This setting determines if the mean or median of those matched peaks should be used as the representative spectrum matching ratio.

## Extras
### MGF Library Prep

The MGF library used to develop zoDIAq lacked protein references. 
We therefore created a function that adds protein references from a given FASTA file. 
While not part of the original publication and therefore not in the GUI, we kept the functionality as made it accessible from the command line.
`zodiaq mgf -m <path to mgf file> -f <path to fasta file>`

Additionally, the `-i` tag can limit peaks to those identified as a fragment of the peptide sequence. This setting was not used in the publication.

## Citations
