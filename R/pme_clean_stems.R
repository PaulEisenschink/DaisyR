#' Clean Stems
#'
#' Cleans up the segmented stems by pme_stem_segm, pme_stem_segm_cell & pme_stem_segm_dbscan
#' @param stems Stem product from pme_stem_segm, pme_stem_segm_cell, or pme_stem_segm_dbscan
#' @param th_nn Maximum of 50cm height intervals which can be empty, before it is no longer considered a stem
#' @param th_area Maximum area of vertically flattened stem points, before it is no longer considered a stem
#' @param th_pts Minimum number of points in the stem, before it is no longer considered a stem
#' @param crs EPSG-code of CRS in which the data is present
#' @return stems_clean List of sf-dataframe & cleaned pointcloud of stems. Separate with st_as_sf(stems_clean[1:5], coords = c('center_x', 'center_y')) & stems_clean[[6]]
#' @export
pme_clean_stems <- function(stems, th_nn = 12, th_area = 1.5, th_pts = 100, crs){
  require(sf)
  require(lidR)
  require(progress)
  
  #th_nn # maxmimum lentgh of 50 cm intercalls which can be empty; lower should remove stricter
  #th_pts # min points to be registered as stem
  #th_area # th of max area in 2d by stem
  
  #create sdf from input pme_stem_segm
  stems_df <- data.frame(stems[1:3])
  stems_sdf <- st_as_sf(stems_df, coords = c('center_x', 'center_y'))
  st_crs(stems_sdf) <- crs
  
  stems_las <- stems[[4]]
  
  #create breaks for nn calc
  i = 1
  hmax <- round(max(stems_las$Z))
  br <- seq(0.5, hmax, by = 0.5)
  
  number_all <- length(stems_sdf$stemID)
  
  pb <- progress_bar$new(total = number_all, format = "  Cleaning Stems [:bar] :percent eta: :eta",)
  
  #get stem basic data: nn, area, len points
  for (sID in stems_sdf$stemID){
    stem_las <- filter_poi(stems_las, stemID == sID)
    len_xyz <- length(stem_las$X)
    stem_xyz <- point_metrics(stem_las, k = 2, xyz = TRUE, func = mean(Z))
    h <- hist(stem_xyz$Z, plot = FALSE, breaks = br)
    stem_nn <- sum(h$counts == 0)
    stems_sdf$nn[i] <- stem_nn
    stems_sdf$area[i] <- st_area(stem_las)
    stems_sdf$len[i] <- len_xyz
    i = i + 1
    pb$tick()
  }
  
  #filter out false trees
  stems_sdf_clean <- stems_sdf[stems_sdf$nn < th_nn, ]
  stems_sdf_clean <- stems_sdf_clean[stems_sdf_clean$area < th_area, ]
  stems_sdf_clean <- stems_sdf_clean[stems_sdf_clean$len > th_pts, ]
  
  #we have to return to df to obtain both df and las
  #stems_df_clean <- data.frame(stems_sdf_clean)
  
  #only take rel stems in las
  las_stems_rel <- filter_poi(stems_las, stemID %in% stems_sdf_clean$stemID)
  
  return(c(stems_sdf_clean, las_stems_rel))
}