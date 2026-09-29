# Potential source areas for atmospheric lead reaching Ny-Ålesund

This repository contains the **code and supplementary material** supporting the study:

> Bazzano, A.; Bertinetti, S.; Ardini, F.; Cappelletti, D.; Grotti, M. *Potential Source Areas for Atmospheric Lead Reaching Ny-Ålesund from 2010 to 2018*. **Atmosphere** 2021, 12, 388.

The work combines **PM10 concentrations, lead isotope ratios, statistical analysis, and atmospheric back-trajectory information** to investigate potential source areas contributing to atmospheric lead measured at Ny-Ålesund, Svalbard, between 2010 and 2018.

The repository is primarily a **research reproducibility archive**: it preserves the R code and supplementary outputs used to reproduce the statistical analyses and the main figures associated with the paper.

## Scientific context

Atmospheric lead at Arctic monitoring sites can originate from sources far from the sampling location.

The analysis therefore combines two complementary types of information:

- **chemical measurements**, including PM10 and lead isotope ratios;
- **atmospheric transport information**, used to investigate the geographical context of the observations.

The statistical analysis explores the distribution and relationships within the measured dataset and supports the interpretation of potential source areas.

The repository does **not** contain the complete back-trajectory analysis workflow. The code, data and results associated with that part of the study are outside the scope of this archive.

## Data and reproducibility

The main PM10 and lead-isotope dataset used in the study is archived separately on **Zenodo**:

**DOI:** 10.5281/zenodo.4484137

The repository therefore separates the reproducibility components into:

1. source code for the statistical analysis;
2. supplementary material and generated results;
3. the archived analytical dataset;
4. the published scientific article.

This makes it possible to distinguish the code preserved here from the external data resources required to reproduce the analysis.

## Repository workflow

The main analysis is contained in:

```text
script.R
```

The script:

1. loads the required R packages;
2. creates the local `dataset/` and `output/` directories;
3. downloads or prepares the required input data;
4. defines functions used throughout the analysis;
5. performs the statistical analyses;
6. reproduces numerical and textual results reported in the manuscript;
7. generates the main figures and saves them as PDF and PNG files.

The generated figures are written to:

```text
output/
```

The workflow is intentionally close to the analysis used for the publication, so that the connection between code, results, figures and manuscript sections remains explicit.

## Requirements

The original analysis was developed and tested with **R 4.0.3**.

The code depends on the following packages:

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

This is a research archive from the original publication rather than a newly structured R package. The dependency list is therefore preserved to support reproduction of the historical analysis.

## Running the analysis

Clone the repository and open `script.R` in R.

The input data should be available in:

```text
dataset/
```

The analysis can then be run from the repository directory.

On successful completion, numerical and textual results are printed to the R session and the reproduced figures are saved in:

```text
output/
```

The exact appearance of figures may depend on the R version and package versions used to run the historical code.

## Publication and citation

If you use the code or supplementary material, please cite the original article:

> Bazzano, A.; Bertinetti, S.; Ardini, F.; Cappelletti, D.; Grotti, M. Potential Source Areas for Atmospheric Lead Reaching Ny-Ålesund from 2010 to 2018. *Atmosphere* 2021, 12, 388. https://doi.org/10.3390/atmos12030388

The dataset used in the analysis is archived separately on Zenodo:

https://doi.org/10.5281/zenodo.4484137

## Supplementary material

Additional results are available in:

```text
supplementary-material.pdf
```

The repository also includes the graphical abstract used to describe the study.

## Scope and limitations

This repository should be read together with the published article.

It provides the source code and supplementary material for the statistical component of the study, but it does not reproduce every element of the broader atmospheric transport analysis. In particular, **back-trajectory data, code and results are not included here**.

The repository is therefore best understood as a reproducibility archive for a published research analysis rather than as a general-purpose workflow for atmospheric source apportionment.

## License

The code is released under the **GNU General Public License**. See `LICENSE.txt` for details.
