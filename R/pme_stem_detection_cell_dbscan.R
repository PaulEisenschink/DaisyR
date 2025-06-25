#' Stem detection algorithm based on DBSCAN & only
#' 
#' Segments stems by cells in height-normalised LiDAR data, based on two step DBSCAN clustering, works singificantly faster then pme_stem_detection_cell
#' 
#' @param las Pointcloud data to be segmented
#' @param minPts Minimum points for the HDBSCAN on initial clustering to be considered as a stem. Higher values disregards more less dense clusters
#' @param cellsize Cell size of initial clustering. Higher cell sized will require considerably more RAM! CAUTION!!!
#' @param hmin Minimum height of point cloud. Points below are disregarded in segmentation
#' @param hmax Maximum height of point cloud. Points above are disregarded in segmentation
#' @param minPts_post Minimum points for the DBSCAN on secondary clustering to be considered as a stem. Higher values disregards more less dense clusters
#' @param eps Epsilon of DBSCAN. Maximum distance between two points to be considered part of the stem. Smaller values segment stricter, but might remove some stem points
#' @export
pme_stem_detection_cell_dbscan <- function(las, minPts = 100, cellsize = 25, hmax = 12, hmin = 0.5, minPts_post = 125, eps = 0.12){
  require(lidR)
  require(dbscan)
  require(stats)
  require(Rdimtools)
  require(progress)
  if (min(las$Z) > 5){
    print('Error:   Point Cloud not normalised!')
  } else {
    cells_done <- 0
    seq_ns <- seq(from = min(las$Y) - 1, to = max(las$Y) + 10, by = cellsize)
    seq_ew <- seq(from = min(las$X) - 1, to = max(las$X) + 10, by = cellsize)
    cells_total <- length(seq_ns) * length(seq_ew)
    pb <- progress_bar$new(total = cells_total, format = "  Segmenting Stems [:bar] :percent eta: :eta",)
    stems_df <- data.frame()
    las_stems <- las[1:5, ]
    dummy_clus <- hdbscan(cbind(las_stems$X, las_stems$Y), minPts = 2)
    las_stems <- add_attribute(las_stems, dummy_clus$cluster, 'cluster')
    for (step_ns in seq_ns){
      for (step_ew in seq_ew){
        cells_done <- cells_done + 1
        #print('Clustering cell...')
        #step_ns <- seq_ns[1]
        #step_ew <- seq_ew[3]
        #print(step_ew)
        #print(step_ns)
        pb$tick()
        las_cell <- filter_poi(las, Y < step_ns + cellsize & Y > step_ns & X < step_ew + cellsize & X > step_ew & Z < hmax & Z > hmin)
        #plot(las_cell)
        if (length(las_cell$X) > 50){
          las_cluster <- dbscan(cbind(las_cell$X, las_cell$Y), minPts = minPts, eps = 0.15)
          las_cell <- add_attribute(las_cell, las_cluster$cluster, 'cluster')
          #plot(las_cell, color = 'cluster')
          #print('Filtering stems...')
          for(clus in unique(las_cell$cluster)){
            clus_single <- filter_poi(las_cell, cluster == clus)
            if (length(clus_single$X) > 10){
              pca_clus <- prcomp(cbind(clus_single$X, clus_single$Y, clus_single$Z))
              pca1 <- abs(pca_clus$rotation[1,1])
              pca2 <- abs(pca_clus$rotation[2,1])
              pca3 <- abs(pca_clus$rotation[3,1])
              #fd_clus <- fd(clus_single$X, clus_single$Y, clus_single$Z)
              z_data_hist <- hist(clus_single$Z, plot = FALSE)
              z_data_sd <- sd(z_data_hist$density)
              z_data_low <- z_data_hist$counts < 5
              z_data_low_sum <- length(z_data_low[z_data_low == TRUE])
              length_stem <- max(clus_single$Z) - min(clus_single$Z)
              if (pca3 > 0.95 & pca1 + pca2 < 0.3 & z_data_low_sum < 5 & z_data_sd < 0.8 & length_stem > 6){ #pca3 > 0.95 & pca1 + pca2 < 0.2 & z_data_low_sum < 5 & z_data_sd < 0.8 & length_stem > 6
                x_cent <- mean(clus_single$X)
                y_cent <- mean(clus_single$Y)
                las_stems <- rbind(las_stems, clus_single)
                stems_df <- rbind(stems_df, c(x_cent, y_cent, z_data_sd, z_data_low_sum, pca3, pca2, pca1, length_stem))
              } else {
                #print(clus)
              }
            }
          }
        }
      }
    }
    print('Finishing up...')
    #drop dummy points
    las_stems <- las_stems[-(1:5), , drop = FALSE]
    #final clustering with dbscan
    nach_clus <- dbscan(cbind(las_stems$X, las_stems$Y), minPts = minPts_post, eps = eps)
    las_stems <- add_attribute(las_stems, nach_clus$cluster, 'stemID')
    las_stems <- filter_poi(las_stems, stemID != 0)
    #get center coords and stemID 
    stems_df <- data.frame()
    for (stamm in unique(las_stems$stemID)){
      stamm_las <- filter_poi(las_stems, stemID == stamm)
      center_x <- mean(stamm_las$X)
      center_y <- mean(stamm_las$Y)
      stamm_data <- c(center_x, center_y, stamm)
      stems_df <- rbind(stems_df, stamm_data)
    }
    colnames(stems_df) <- c('center_x', 'center_y', 'stemID')
    return(c(stems_df, las_stems))
  }
}



