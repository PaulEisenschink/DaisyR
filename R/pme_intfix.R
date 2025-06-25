#' GS-100C+ Intensity with flightlines
#' 
#' Corrects the flight lines effect from raw, stitched GS-100C+ LiDAR data
#' 
#' @param input_path Path to LAS file
#' @param drone_path Path to GNSS flight file 
#' @param crs EPSG code of the LiDAR data
#' @return NULL Writes the corrected LAS File into the input directory 
#' @export
pme_intfix <- function(input_path, drone_path, crs){
  file_no_ext <- gsub('.{4}$', '', input_path)
  print("reading input file...")
  las_cat = read.las(input_path)
  print("fixing point source...")
  ##PointSourceID entfernen
  las_cat = las_cat %>% select(-PointSourceID)
  
  ##Neuen Vektor mit PointSourceID erstellen
  PointSourceID = rep(2, times = length(las_cat$X))
  ##PointSourceID muss integer sein, deshalb as.integer verwenden
  PointSourceID = as.integer(PointSourceID)
  
  ##Neuen Vektor als Spalte anh?ngen
  las_cat = cbind(las_cat, PointSourceID)
  
  ##Header erstellen
  las_head = header_create(las_cat)
  
  ##las File erstellen
  write.las(paste0(file_no_ext, "_temp.las"), las_head, las_cat)
  
  gc()
  
  ######################################
  #read las file
  las_fix = readLAS(paste0(file_no_ext, "_temp.las"))
  print("reading trajectory data...")
  #read GNSS data from LiDAR processing
  drone_lonlat = read.csv(drone_path, sep = '', header = FALSE)
  
  las_fix@data$gpstime[1]
  
  las_fix@data$gpstime[length(las_fix@data$gpstime)]
  
  sensor_1 = drone_lonlat[drone_lonlat$V2 >= las_fix@data$gpstime[1],]
  
  sensor_1 = sensor_1[,-1]
  names(sensor_1)[names(sensor_1) == "V3"] <- "X"
  names(sensor_1)[names(sensor_1) == "V4"] <- "Y"
  names(sensor_1)[names(sensor_1) == "V5"] <- "Z"
  names(sensor_1)[names(sensor_1) == "V2"] <- "gpstime"
  sensor_1 = sensor_1[,0:4]
  
  #Richtige Pointsourceid einfügen
  PointSourceID = rep(7, times = nrow(sensor_1))
  sensor_1$PointSourceID = PointSourceID
  XYZ = sensor_1[,c(2,3,4)]
  sensor_2 = SpatialPointsDataFrame(data = sensor_1, coords = XYZ)
  proj4string(sensor_2) <- CRS("+proj=longlat +datum=WGS84")
  sensor_2 <- spTransform(sensor_2, crs(paste0("EPSG:", crs)))
  
  print("normalizing intensities...")
  #sensor = track_sensor(las1, Roussel2020(interval = 6, pmin = 1))
  las_fix <- normalize_intensity(las_fix, range_correction(sensor_2, Rs = 50))
  print("writing fixed file...")
  #LAS speichern 
  writeLAS(las_fix, paste0(file_no_ext, "_fix.laz"))
  print("cleaning up...")
  file.remove(paste0(file_no_ext, "_temp.las"))
  gc()
  print("Done!")
}