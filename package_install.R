# Installation of required R packages --------------------------------

# required packages
packages <- c("sf",
              "raster",
              "rgdal",
              "cowplot",
              "geosphere",
              "pracma",
              "progress",
              "plot3D",
              "GWmodel",
              "psych",
              "roahd",
              "ggplot2",
              "reshape2")

# new packages
new_packages  <- packages[!(packages %in% installed.packages()[,"Package"])] 

# install required packages
if(length(new_packages))
{
  install.packages(new_packages)
}