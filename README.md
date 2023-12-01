# zoDIAq
## Introduction

*Uses an MIT License* 

zoDIAq (derived from "csodiaq" acronym, Cosine Similarity Optimization for DIA qualitative and quantitative analysis) is a software tool 
for the analysis of direct infusion shotgun proteome analysis (DISPA) data and data independent acquisition (DIA) data.

[Click here for the paper published in *Analytical Chemistry*.](https://doi.org/10.1021/acs.analchem.1c02021)

![Paper Title Header](https://github.com/xomicsdatascience/zoDIAq/blob/main/img/zodiaq-paper.png)

## Instructions
### Installation

All installation instructions are done from the command line. Download and use of zoDIAq requires access to the `git` and `pip` packages. For system-specific instructions on using the command line or installing these packages, see the [wiki page.](https://github.com/xomicsdatascience/zoDIAq/wiki)

* [Windows install instructions click here](https://github.com/xomicsdatascience/zoDIAq/wiki/Install-Instructions-for-Windows)
* [Mac install instructions click here](https://github.com/xomicsdatascience/zoDIAq/wiki/Install-Instructions-for-Mac)

### Accessing the GUI

From the command line enter `zodiaq gui` to start the GUI. All instructions in this README file pertain to use of the GUI. 
Instructions for use from the command line can be accessed by entering `zodiaq -h`.

### Peptide Identification

![Peptide Identification GUI Picture](https://github.com/xomicsdatascience/zoDIAq/assets/11773171/6356880f-80df-4baa-83c4-ac2675108ff0)

Files: 
1. `DIA Data File`: This is a required field. Choose at least one data file from a DISPA or DIA data run.
2. `Library File`: This is a required field. Choose a reference library file that is in TraML (.tsv or .csv) or MGF format. 
zoDIAq treats .tsv files as tab-deliminated and .csv as comma-deliminated. Some library formats do not follow this standard (ie, .csv files that are tab-deliminated) and will require adjustment by the end user.
3. `Outfile Directory`: This is a required field. Choose a folder that the output files should go into.

Settings:

Note that it is recommended that the default value is used for each of these settings.

1. `Pre-correction m/z difference tolerance for peak matches (in ppm)`: Default value is 30. 
If you're not sure what would make an ideal setting, go ahead and leave this blank, but check the histogram box below. 
If the resulting histogram doesn't have a normal distribution (a peak), consider using a wider tolerance.
2. `Correction`: Default is checked. This value resets the PPM tolerance range based on an initial scan. It reduces the liklihood of identifying false positives.
3. `Corrective Standard Deviations`: Only available if the correction box is checked. Default is to customize the corrected tolerance to the distribution of the histogram, excluding noise around the peak. 
Entering a value here instead sets the tolerance to a standard deviation of the distribution.
4. `Create Histogram`: Only available if the correction box is checked. Default is checked. 
Generates a histogram to visualize the corrective PPM tolerance used in the corrected analysis.

### Scoring and FDR Filtering

![Scoring and FDR Filtering](https://github.com/xomicsdatascience/zoDIAq/assets/11773171/c0c36ffa-676c-4a7c-b30b-0fb66517d1d6)

Files: 
1. `zoDIAq Scoring Output Directory`: This is a required field. Choose the directory created by the `Peptide Identification` output of zoDIAq.

Settings:

None.

### Targeted Peptide Reanalysis

![Screenshot 2023-12-01 at 10 25 53 AM](https://github.com/xomicsdatascience/zoDIAq/assets/11773171/7495fcfc-65ae-4197-b092-1624c88bf18e)

Files: 
1. zoDIAq Identification Output Directory: This is a required field. Choose the directory created by the `Scoring and FDR Filtering` output of zoDIAq.

Settings:

1. `Maximum Number of peptides per protein`: By default, no protein analysis is completed. If a number is provided, the program will provide a mass spectrometer file targetting the given number of peptides for each protein.
2. `Proximity m/z values should be to a bin value`: To avoid increased cost for targetting every peptide individually, specific m/z bins are targeted for overlapping peptides. This setting indicates how close a peptide's m/z value should be to the "center" of a bin.
3. `Include Heavy Isotopes for SILAC protocol`: For targeted reanalysis of peptides using a SILAC protocol, the heavy isotopes for chosen peptides will be included as well.

## Citations
