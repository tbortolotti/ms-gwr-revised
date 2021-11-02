# ms-gwr-revised

## Multi-Source Geographically Weighted Regression for Regionalized Ground-Motion Models - A Code Review

This is a commented review of the GitHub repository of Luca Caramenti, which contains the supplementary material for Multi-Source Geographically Weighted Regression (MS GWR)

[lucaramenti/ms-gwr](https://github.com/lucaramenti/ms-gwr)

The purpose of a revised version of the code is to provide a user-friendly tool which is easily accessible by users who are new to the subject of MS GWR.

## Structure of the repository
The repository is composed as follows:
* `simulated_data_example.R`: File containing the simulation study conducted by Luca Caramenti.
* `main.R`: R script containing the MS GWR analysis on the Italian seismic data for the PGA.
* `data_dir`: Folder containing the data of the case study. In particular, `confini_ut33.shp` is the shape file for the Italian territory for which the analysis is performed, and `italian_data_pga.RData` contains the full dataset for the PGA measurements.
* `functions`: Folder containing all functions used in the analysis, both for the simulation study and the case study. Each function in the folder contains a brief description of its usage and of its input and return parameters.

### Notebook for the application to Italian seismic data
[Report of MS GWR application to Italian seismic data](https://htmlpreview.github.io/?https://github.com/lucaramenti/ms-gwr/blob/main/msgwr_seismological_data_notebook.nb.html)

This is a useful notebook reporting the entire procedure applied by Luca Caramenti to the dataset taken from Engineering Strong Motion database
> Luzi, L., G. Lanzano, C. Felicetta, M.C. D'Amico, E. Russo, S. Sgobba, F. Pacor and ORFEUS Working Group 5. Engineering strong motion database (ESM), 2020. URL: [https://esm-db.eu](https://esm-db.eu).

This procedure is 

### Installation

Run file `package_install.R` to have an automatic installation of the required R packages.

## Author
Teresa Bortolotti
