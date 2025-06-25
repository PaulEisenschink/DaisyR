#pipeline: pme_stem_segmentation_dbscan > pme_clean_stems > pme_dbh_estimation


pme_dbh_estimation <- function(stems_clean, z_mid, q, crs){
  stems_las <- stems_clean[[6]]
  stems_df <- data.frame(stems_clean[1:5])
  stems_sdf <- st_as_sf(stems_df)
  st_crs(stems_sdf) <- crs
  
  dbh_df <- data.frame()
  for (ID in stems_sdf$stemID){
    stem_single <- filter_poi(stems_las, stemID == ID)
    stem_single <- filter_poi(stem_single, Z > (z_mid - 1) & Z < (z_mid + 1))
    if (length(stem_single$X > 20)){
      pca <- prcomp(cbind(stem_single$X, stem_single$Y, stem_single$Z))$rotation
      pca_th <- abs(pca[2,1]) + abs(pca[1,1])
      if (pca_th < 0.7){
        stem_single_xy <- point_metrics(stem_single, xyz = TRUE, func = mean(Z), k = 2)
        stem_single_xy_d <- dist(cbind(stem_single_xy$X, stem_single_xy$Y))
        x_coords <- mean(stem_single_xy$X)
        y_coords <- mean(stem_single_xy$Y)
        tryCatch(
          expr = {dbh_df <- rbind(dbh_df, c(ID, x_coords, y_coords, quantile(stem_single_xy_d, q)))},
          error = function(e){
            dbh_df <- rbind(dbh_df, c(ID, NA, NA, NA))
          }
        )
      }
    } else {
      dbh_df <- rbind(dbh_df, c(ID, NA, NA, NA))
    }
  }
  
  colnames(dbh_df) <- c('ID', 'X', 'Y', 'dbh_algo')
  dbh_df <- dbh_df[!is.na(dbh_df$X), ]
  dbh_sdf <- st_as_sf(dbh_df, coords = c('X', 'Y'), crs = crs)
  
  return(dbh_sdf)
}
