# ms-gwr-revised

## Multi-Source Geographically Weighted Regression for Regionalized Ground-Motion Models - A Code Review

This is a commented review of the GitHub repository of Luca Caramenti, which contains the supplementary material for Multi-Source Geographically Weighted Regression (MS GWR)

[lucaramenti/ms-gwr](https://github.com/lucaramenti/ms-gwr)

The purpose of a revised version of the code is to provide a user-friendly tool which is easily accessible by users who are new to the subject of MS GWR.

## Structure of the repository
The repository is composed as follows:
* `simulated_data_example.R`: File containing the simulation study conducted by Luca Caramenti.
* `main.R`: R script containing the MS GWR analysis on the Italian seismic data for the PGA, with dependence on site and event coordinates.
* `main_midpoint.R`: R script containing the MS GWR analysis on the Italian seismic data for the PGA, with dependence on site coordinates and midpoint of site and event coordinates.
* * `data_preparation.R` is the code used to perform data preprocessing. It loads files `data_dir/confini_ut33.shp` and `data_dir/italian_data_pga.RData` and prepares data `data_dir/utm_coordinates.RData` and `data_dir/regressors.RData` that are loaded in `main.R` for the analysis.
* `data_dir`: Folder containing the data of the case study. In particular, `confini_ut33.shp` is the shape file for the Italian territory for which the analysis is performed, and `italian_data_pga.RData` contains the full dataset for the PGA measurements. `utm_coordinates.RData` and `regressors.RData`contain preprocessed data, that are loaded to run the analysis.
* `functions`: Folder containing all functions used in the analysis, both for the simulation study and the case study. Each function in the folder contains a brief description of its usage and of its input and return parameters.

## Analysis of Italian seismic data
The revised version of the code focuses mainly on the application of MS GWR for the analysis of Italian seismic data. File `main.R` builds the optimal model for the PGA (Peak Ground Acceleration), conditionally on seismic parameters of magnitude (mag), shear-wave velocity (vs30), Joyner_Boore distance (JB_complete) and style-of-faulting (fm_type_code). Included in the analysis are the longitude and the latitude of the seismic events and of the sites of registration of the seismic waves.

The course-of-action of the analysis may be sketched as follows:

0. Preprocessing: Removal of undesired sites and events and dataset building.
1. Selection of the optimal bandwidth for event- and site-related gaussian kernels.
2. Permutation tests for assessing the non-stationarity of regression coefficients.
3. Permutation tests for the significance of constant regression coefficients.
4. GCV comparison between SEC and ESC algorithms, to see which algorithm performs better with the data.
5. Computation of R<sup>2</sup><sub>adj</sub>.
6. Full calibration of the model and computation of the regression coefficients over a spatial grid of interest.

The procedure closely follows what is reported in

[Report of MS GWR application to Italian seismic data](https://htmlpreview.github.io/?https://github.com/lucaramenti/ms-gwr/blob/main/msgwr_seismological_data_notebook.nb.html)

, which is a useful notebook reporting the analysis of Luca Caramenti on data taken from Engineering Strong Motion database
> Luzi, L., G. Lanzano, C. Felicetta, M.C. D'Amico, E. Russo, S. Sgobba, F. Pacor and ORFEUS Working Group 5. Engineering strong motion database (ESM), 2020. URL: [https://esm-db.eu](https://esm-db.eu).

## Installation

Run file `package_install.R` to have an automatic installation of the required R packages.

### Author
Teresa Bortolotti
