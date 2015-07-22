#' Download Requried Packages
#' 
#' Automatically downloads all the required packages for full analysis
#' @export
download_packages <- function() {
  
  try(install.packages("lme4", dependencies=T), F)
  
  try(install.packages("multcomp", dependencies=T), F)
  
  try(install.packages("seqinr", dependencies=T), F)
  
  try(install.packages("qpcR", dependencies=T), F)
  
}



