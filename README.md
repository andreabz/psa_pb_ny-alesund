# Supplementary material for Potential source areas for atmospheric lead reaching Ny-Ålesund from 2010 to 2018.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4554221.svg)](https://doi.org/10.5281/zenodo.4554221)

### Authors:
Andrea Bazzano<sup>1,</sup>\*,
Stefano Bertinetti<sup>1</sup>,
Francisco Ardini<sup>1</sup>,
David Cappelletti<sup>2</sup> and 
Marco Grotti<sup>1</sup>.

<sup>1</sup> Department of Chemistry and Industrial Chemistry, University of Genoa, via Dodecaneso 31, 16146, Genoa, Italy

<sup>2</sup> Department of Chemistry, Biology and Biotechnologies, University of Perugia, Via Elce Di Sotto 8, 06123, Perugia, Italy

\* corresponding author: andrea dot bazzano at edu dot unige dot it, orcid id: https://orcid.org/0000-0002-9353-3919

This repository contains source code and additional supplementary materials from our manuscript, "Supplementary material for Potential source areas for atmospheric lead reaching Ny-Ålesund from 2010 to 2018". Additional results can be found within supplementary-material.pdf. The main dataset with PM<sub>10</sub> and lead isotope ratio measured values has been archived on Zenodo for reproducibility (http://doi.org/10.5281/zenodo.4484137).
Data, code and results for back-trajectory analysis are not included in this repository.

The following instructions provide details on how to run the source code underlying the analysis, including replication of the main figures and results.

## Requirements
The code has been tested with R version 4.0.3, "Bunny-Wunnies Freak Out" The following R packages and their dependencies must be installed before the code will run successfully:

- `data.table`
- `dplyr`
- `lubridate`
- `summarytools`
- `fitdistrplus`
- `dunn.test`
- `mclust`
- `QuantPsyc`
- `energy`
- `MASS`
- `ggplot2`
- `ggpubr`
- `ggrepel`
- `ggforce`
- `patchwork`
- `scales`

## Instructions

Before running the code, make sure the required R packages have been installed.  Set the R working directory to the location of this README file. Input data need to be downloaded and saved in the `./dataset/` subdirectory of the R working directory, whereas figures generated running the code will be saved in the `./output/` subdirectory of the R working directory.

Running the entire script will require few minutes on most computers.

### Step One: 

- Loads required packages into R.

### Step Two: 

- Creates `./dataset/` and `./output/` subdirectories, which will hold the underlying data sources and analysis output, respectively.

### Step Three:

- Downloads the raw data sources used in the analysis. These data are publicly available in the repository and the main dataset has been archived on Zenodo for reproducibility (doi = http://doi.org/10.5281/zenodo.4484137). Input data require less than 50 kB of space.

### Step Four: 

- Open the script `script.R` in R and run the code. The script start defining some functions used in the following analysis. Textual and numerical results are presented citing the sections of the submitted manuscript. Tables and Figures are reproduced at the end of the script. Only figures are saved in `./output/` subdirectory.

## Output

Upon successful completion of `script.R`, numerical and textual results are visualized on screen and figures are saved as PDFs and PNGs in the `./output/` folder. Example format includes `./output/figure1.png`, etc. These figures will look very similar, if not identical, to those found in the manuscript.

## License
The entire code is available under the GNU General Public license. See LICENSE.txt for more information
