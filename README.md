<img width="1767" height="917" alt="image" src="https://github.com/user-attachments/assets/233f58dd-1c3a-4f0a-818e-5387e3d36b37" />














# TAGEM-GDD-Dynamics

## Long-Term Growing Degree Day Dynamics Across TAGEM Research Locations in Türkiye (1981–2025)

This repository contains the data inputs and R scripts used to characterize long-term changes in growing degree day (GDD) dynamics across 30 research locations of the General Directorate of Agricultural Research and Policies (TAGEM) in Türkiye during the 1981–2025 period.

Daily meteorological data are retrieved from the NASA POWER database using the geographic coordinates of the research locations. Growing degree days are calculated using a base temperature of 10 °C.

The analytical workflow combines longitudinal functional data analysis, temporal trend analysis, functional principal component analysis, k-means clustering, and spatial visualization to characterize long-term changes in thermal accumulation across TAGEM research environments.

---

## Repository Contents

### `Institutes_with_coords.xlsx`

Contains the identifiers and geographic coordinates (latitude and longitude) of the TAGEM research locations included in the study.

The coordinates are used as spatial inputs for retrieving location-specific daily meteorological data from the NASA POWER database.

### `Download_NASA_POWER_1981_2025_Save_RDS.R`

Retrieves daily meteorological data from NASA POWER for the selected TAGEM research locations over the 1981–2025 period using the `envirotypeR` package.

The downloaded meteorological data are stored in RDS format and used as input for the subsequent GDD analyses.

### `LFPCA_FPCA_Kmeans.R`

Contains the main analytical workflow for characterizing long-term GDD dynamics across the research locations.

The analysis consists of two complementary functional frameworks.

---

## Analytical Framework

### 1. Longitudinal Functional Analysis of Daily GDD

Daily GDD is calculated using a base temperature of 10 °C:

**GDD = max[((TMAX + TMIN) / 2) − 10, 0]**

Annual daily GDD trajectories are constructed for each research location and year.

Longitudinal functional principal component analysis (LFPCA) is then used to characterize variation in these trajectories across both the seasonal (day-of-year) and longitudinal (year) dimensions.

The main steps include:

- construction of environment × year daily GDD trajectories;
- estimation of the longitudinal GDD mean surface;
- marginal functional principal component analysis;
- estimation of location-specific longitudinal score trajectories;
- secondary principal component analysis of longitudinal features;
- dimensionality reduction based on cumulative explained variance;
- k-means clustering of research locations;
- evaluation of cluster structure; and
- spatial visualization of the resulting long-term thermal clusters.

### 2. Daily Long-Term GDD Trend Analysis

A complementary analysis is conducted to quantify how daily GDD has changed over the 45-year study period.

For each research location and each day of the year, a linear temporal model is fitted across years:

**GDD = Intercept + Slope × (Year − 1981)**

The estimated slope therefore represents the long-term annual rate of GDD change for a particular research location and day of the year.

Together, the 365 daily slopes form a complete annual GDD-change trajectory for each research location.

The analysis includes:

- estimation of daily long-term GDD trend slopes;
- calculation of model statistics;
- characterization of seasonal patterns in GDD change;
- monthly aggregation of daily trend slopes;
- calculation of monthly integrated slope statistics using trapezoidal area under the curve (AUC); and
- comparison and ranking of research locations according to monthly GDD-change patterns.

### 3. Functional Analysis of Daily GDD Trend Slopes

Functional principal component analysis (FPCA) is applied to the 365-day long-term GDD trend-slope trajectories.

This framework summarizes differences among research locations in both the magnitude and seasonal timing of long-term GDD change.

The main steps include:

- FPCA of daily GDD trend-slope curves;
- selection of functional principal components based on cumulative variance explained;
- extraction of location-specific FPCA scores;
- visualization of functional eigenfunctions;
- k-means clustering based on retained FPCA scores;
- evaluation of alternative clustering solutions using silhouette statistics; and
- spatial mapping of the resulting GDD-change clusters across Türkiye.

---

## General Workflow

`TAGEM research location coordinates`

↓

`NASA POWER daily meteorological data (1981–2025)`

↓

`Daily GDD calculation (Tbase = 10 °C)`

↓

**Longitudinal framework**

`Annual daily GDD trajectories → LFPCA → Secondary PCA → k-means clustering`

↓

**Temporal-change framework**

`Daily GDD trends → 365-day slope trajectories → FPCA → k-means clustering`

↓

`Monthly integrated GDD trend analysis`

↓

`Spatial characterization of long-term thermal dynamics across Türkiye`

---

## Main R Packages

The workflow relies primarily on the following R packages:

- **`envirotypeR`** — retrieval of daily meteorological data from the NASA POWER database for the selected research locations.
- **`refund`** — longitudinal functional principal component analysis (LFPCA) of daily GDD trajectories.
- **`fdapace`** — functional principal component analysis (FPCA) of daily long-term GDD trend-slope trajectories.
- **`sf`** — spatial processing and mapping of TAGEM research locations and derived clusters.

Additional R packages for data manipulation and visualization are specified directly within the corresponding scripts.

---

## Data Source

Daily meteorological data are obtained from the NASA Prediction Of Worldwide Energy Resources (NASA POWER) database.

The latitude and longitude of the TAGEM research locations required for data retrieval are provided in `Institutes_with_coords.xlsx`.

Because the meteorological data can be retrieved directly from NASA POWER using the provided coordinates and R script, the complete downloaded weather dataset does not need to be stored in this repository.

---

## Reproducibility

The scripts should be executed in the following order:

1. `Download_NASA_POWER_1981_2025_Save_RDS.R`
2. `LFPCA_FPCA_Kmeans.R`

The first script retrieves the daily meteorological data from NASA POWER and stores them in RDS format.

The second script calculates daily GDD and performs the longitudinal functional, temporal trend, FPCA, clustering, and spatial analyses.

File paths may need to be adjusted according to the user's local working directory before running the scripts.

---

## Software

All data processing and statistical analyses were conducted in **R**.

Package requirements and analytical settings are provided within the corresponding R scripts to facilitate reproducibility.

---

## Citation

If you use the data-processing workflow, analytical methods, or R code provided in this repository, please cite the associated publication.

**Citation information will be added following publication.**

---

## Contact

**Alper Adak**  alp.adk.07@gmail.com
Batı Akdeniz Agricultural Research Institute (BATEM), 
General Directorate of Agricultural Research and Policies (TAGEM), 
T.C. Ministry of Agriculture and Forestry, 
Antalya, Türkiye (BATEM)
Ankara, Türkiye (TAGEM)
