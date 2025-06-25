#' Integration of Segmented Stem Points
#' 
#' Integrates the segmented stem points from pme_stem_detection_cell or pme_stem_detection_cell_dbscan into the original pointcloud
#' @param las Pointcloud to be reintergrated with the semgented stem points
#' @param las_attr Pointcloud of segmented stems from pme_stem_detection_cell, pme_stem_detection_cell_dbscan, or pme_clean_stems
#' @return Returns Pointcloud with stemID as new Variable: las$stemID

pme_add_stemID <- function(las, las_attr){
  require(lidR)
  las <- add_lasattribute(las, rep(0, length(las$X)), 'stemID', 'stemID')
  las$stemID[las$gpstime %in% las_attr$gpstime] <- las_attr$stemID
  return(las)
}
