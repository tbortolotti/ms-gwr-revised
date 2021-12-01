# ms-gwr-revised

## Multi-Source Geographically Weighted Regression for Regionalized Ground-Motion Models - A Code Review

This is a commented review of the GitHub repository of Luca Caramenti, which contains the supplementary code material for Multi-Source Geographically Weighted Regression (MS GWR)

[lucaramenti/ms-gwr](https://github.com/lucaramenti/ms-gwr)

The purpose of a revised version of the code is to provide a user-friendly tool which is easily accessible by users who are new to the subject of MS GWR and interested in exploiting its properties, both for simulated and seismological applications.

## The Multi-Source Geographically Weighted Regression model
Geographically Weighted Regression (GWR) is a family of satistical methods aimed to estimate a regionalized linear model that is characterized by spatially varying regression coefficients. Motivated by the inherent spatial dependencies of the regression coefficients appearing in the formulation of ground-motion models, GWR models are extended in

> Caramenti, L., A. Menafoglio, S. Sgobba, G. Lanzano. Multi-Source Geographically Weighted Regression for Regionalized Ground-Motion Models, 2020.

to allow for the precence of two-sources of spatial non-stationarity. In particular, event- and station-related regression coefficients enter the regression model as dependent on event and station coordinates, respectively. 
Peculiarity of a GWR model lies in the estimation of regression coefficients that are specific of a point in space, so that each geographical location is coupled to its own vector of coefficients estimates. The introduction of a double dependency on space follows the same rationale, and requires the formulation and implementation of an algorithm that appropriately accounts for the estimation procedure.

Consider a general MS GWR model that allow for the presence of constant, station-dependent and event-dependent coefficients. Then the algorithm, extensively discussed and commented in the work of Caramenti et al. (2020), works in cascade as syntesized in these five steps:
0. Computation of three large matrices H<sub>E</sub>, H<sub>S</sub> and B.
1. Estimation of constant coefficients Beta<sub>C</sub>.
2. Evaluation of the residuals with respect to the constant covariates: res<sub>C</sub> = Y - X<sub>C</sub> Beta<sub>C</sub>.
3. Estimation of station-dependent coefficients Beta<sub>S</sub> using res<sub>1</sub>.
4. Evaluation of the residuals with respect to constant and station-related covariates: res<sub>SC</sub> = res<sub>C</sub> - X<sub>S</sub> Beta<sub>S</sub>.
5. Estimation of event-dependent coefficients Beta<sub>E</sub> using res<sub>SC</sub>.

It is worth noticing that the order of estimation of the coefficients is totally arbitrary. Here a ESC (Event-Station-Constant) procedure is presented, but any permutation of the calibration sequence is possible and leads to equally valid results, that do not critically differ one from the other. The practicioner should test multiple permutations of the algorithm strategy on data under analysis and pick which performs best.

## Structure of the repository
The repository is composed of:
* `simulated_data_example.R`: File containing the simulation study conducted by Luca Caramenti.
* `main_event.R`: R script containing the MS GWR analysis on the Italian seismic data for the PGA, with dependence on station and event coordinates.
* `main_midpoint.R`: R script containing the MS GWR analysis on the Italian seismic data for the PGA, with dependence on station coordinates and station-event midpoint coordinates.
* `data_preparation.R` is the code used to perform data preprocessing. It loads files `data_dir/confini_ut33.shp` and `data_dir/italian_data_pga.RData` and prepares data `data_dir/utm_coordinates.RData` and `data_dir/regressors.RData` that are loaded in `main_event.R` and `main_midpoint.R` for the analysis.
* `data_dir`: Folder containing the data of the case study. In particular, `confini_ut33.shp` is the shape file for the Italian territory for which the analysis is performed, and `italian_data_pga.RData` contains the full dataset for the PGA measurements. `utm_coordinates.RData` and `regressors.RData`contain preprocessed data, that are loaded to run the analysis.
* `functions`: Folder containing all functions used in the analysis, both for the simulation study and the case study. Each function in the folder contains a brief description of its usage and of its input and return parameters. It is relevant to point out that subfolder `simulation` contains the functions used for the simulation study. They perform the same exact analysis than those used for the seismic case study, but are not optimized for the analysis of large datasets. Nonetheless, they are more readable and understandable, so that a practitioner should start from the simulation study and the simulation functions to better understand how they work. This is also the ratio behind their separate presence in this repository.

## Simulation study

The course-of-action of the analysis, that sticks to what is clearly illustrated in

[Report of MS GWR application to simulated data](https://htmlpreview.github.io/?https://github.com/lucaramenti/ms-gwr/blob/main/simulated_data_example_notebook.nb.html)

may be sketched as follows:

0. Preprocessing: Removal of undesired stations and events and dataset building.
1. Selection of the optimal bandwidth for event- and station-related gaussian kernels.
2. Permutation tests for assessing the non-stationarity of regression coefficients.
3. Permutation tests for the significance of constant regression coefficients.
4. GCV comparison between SEC and ESC algorithms, to see which algorithm performs best with the data.
5. Computation of R<sup>2</sup><sub>adj</sub>.
6. Full calibration of the model and computation of the regression coefficients over a spatial grid of interest.

The simulation study explores the performances of MS GWR with particular regard to estimation and prediction accuracy when changing (i) the order of estimation (ESC, SEC, CES, ..) and (ii) the bandwidth of the spatial kernels involved in the GWR estimates. Stationary, station-dependent and event-dependent coefficients are generated and then the procedure displayed above is deployed to obtain their MS GWR estimates.

![alt text](https://github.com/tbortolotti/ms-gwr-revised/blob/main/plots/simulation/coefficient_estimates.png)

## Analysis of Italian seismic data
The revised version of the code focuses mainly on the application of MS GWR for the analysis of Italian seismic data. Files `main_event.R` and `main_midpoint.R` calibrate the optimal model for the PGA (Peak Ground Acceleration), conditionally on seismic parameters of magnitude (mag), shear-wave velocity (vs30), Joyner_Boore distance (JB_complete) and style-of-faulting (fm_type_code). Included in the analysis are the longitude and the latitude of the seismic events, of the stations of registration of the seismic waves, and of the midpoint coordinates of event and station.

Event coordinates                                                                    |  Station coordinates        | Station-Event midpoint coordinates
:-----------------------------------------------------------------------------------:|:-------------------------:  | :-------------------------:
![](https://github.com/tbortolotti/ms-gwr-revised/blob/main/plots/event_coords.png)  |  ![](https://github.com/tbortolotti/ms-gwr-revised/blob/main/plots/site_coords.png)  | ![](https://github.com/tbortolotti/ms-gwr-revised/blob/main/plots/mid_coords.png)

The procedure closely follows what is reported in

[Report of MS GWR application to Italian seismic data](https://htmlpreview.github.io/?https://github.com/lucaramenti/ms-gwr/blob/main/msgwr_seismological_data_notebook.nb.html)

, which is a useful notebook reporting the analysis of Luca Caramenti on data taken from Engineering Strong Motion database
> Luzi, L., G. Lanzano, C. Felicetta, M.C. D'Amico, E. Russo, S. Sgobba, F. Pacor and ORFEUS Working Group 5. Engineering strong motion database (ESM), 2020. URL: [https://esm-db.eu](https://esm-db.eu).

The step forward reported in this repository consists in the optimization of the original code to better handle the analysis on very large datasets. A parallel computation of large matrices is provided, that highly reduces the overall compuational burden.

## Installation

Run file `package_install.R` to have an automatic installation of the required R packages.

### Author
Teresa Bortolotti
