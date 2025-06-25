

pme_understory_detection <- function(las, res = 25, bp_pre = 20, min_uth = 5){
  require(lidR)
  require(terra)
  
  print('Creating Understory Raster...')
  las <- filter_poi(las, stemID == 0 & Z < max(Z)/2 & Z > 0.5 & Classification == 1)
  prelim_bp <- bp_pre
  las <- filter_poi(las, Z < prelim_bp)
  bp_raster <- pixel_metrics(las, ~which.min(hist(Z, breaks = seq(min(las$Z), max(las$Z) + 1, by = 1), plot = FALSE)$counts), res)
  
  bp_raster[values(bp_raster) > 15] <- 0
  bp_raster[values(bp_raster) < 4] <- 5
  
  print('Clipping LAS...')
  i <- 1
  cells_total <- length(values(bp_raster))
  for (cells in values(bp_raster)){
    cell_center_point <- xyFromCell(bp_raster, i)
    cell_lid <- filter_poi(las, X < cell_center_point[1] + res/2 & X > cell_center_point[1] - res/2 & Y < cell_center_point[2] + res/2 & Y > cell_center_point[2] - res/2)
    cell_lid_under <- filter_poi(cell_lid, Z < cells)
    progr <- (i / cells_total) * 100
    print_progress(progr)
    if (i == 1){
      bottom_post <- cell_lid_under
    } else {
      bottom_post <- rbind(bottom_post, cell_lid_under)
    }
    i <- i + 1
  }
  return(bottom_post)
}

