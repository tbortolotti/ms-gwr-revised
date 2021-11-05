#'
#' Installation of the required R packages
#'


# required packages
packages <- c("sf",
              "raster",
              "rgdal",
              "GWmodel",
              "cowplot",
              "geosphere",
              "psych",
              "pracma",
              "reshape2",
              "plot3D",
              "progress",
              "roahd",
              "ggplot2",
              "snowfall")

# new packages
new_packages  <- packages[!(packages %in% installed.packages()[,"Package"])] 

# install required packages
if(length(new_packages))
{
  install.packages(new_packages)
}