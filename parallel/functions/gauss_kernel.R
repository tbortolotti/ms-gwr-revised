#'
#' Gaussian Kernel
#'
#' This function builds a vector of weights as given by a gaussian kernel with bandwidth h
#'
#' @param d:          vector of distances (one column of a distance matrix)
#' @param h:          bandwidth of the kernel
#' 
#' @return wgts:      vector of weights
#' 

gauss_kernel = function(d, h){
  wgts = exp(-0.5*(d/h)^2)
  return(wgts)
}