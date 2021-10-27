# ms-gwr-reviewed

## Multi-Source Geographically Weighted Regression for Regionalized Ground-Motion Models - A Code Review

This is a commented review of the GitHub repository of Luca Caramenti, which contains the supplementary material for Multi-Source Geographically Weighted Regression (MS GWR)

[lucaramenti/ms-gwr](https://github.com/lucaramenti/ms-gwr)

The purpose of a revised version of the code is to provide a user-friendly tool which is easily accessible by users who are new to the subject of MS GWR.

## Structure of the repository
The repository is composed as follows:
* `functions`: Folder containing all functions used in the analysis, both for the simulation study and the case study (se voglio scrivere in grassetto **è così**)
* `simulated_data_example.R`: R-script with all the steps for running the Blocked Gibbs sampler
* `msgwr_seismological_data_notebook.nb.html`: report of the case study analysis developed by Luca Caramenti

### Notebook for the application to Italian seismic data
[Report of MS GWR application to Italian seismic data](https://htmlpreview.github.io/?https://github.com/lucaramenti/ms-gwr/blob/main/msgwr_seismological_data_notebook.nb.html)

This is a useful notebook reporting the entire procedure applied by Luca Caramenti to the dataset taken from Engineering Strong Motion database
> Luzi, L., G. Lanzano, C. Felicetta, M.C. D'Amico, E. Russo, S. Sgobba, F. Pacor and ORFEUS Working Group 5. Engineering strong motion database (ESM), 2020. URL: [https://esm-db.eu](https://esm-db.eu).

### Installation

Run file `package_install.R` to have an automatic installation of the required R packages.

## Author
Teresa Bortolotti
