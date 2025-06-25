pme_tree_curvature <- function(stems, chm = chm, crs = 32633){
  require(dplyr)
  
  #create stems sdf form stems
  stems_df <- data.frame(stems[1:5])
  stems_sdf_clean <- st_as_sf(stems_df)
  st_crs(stems_sdf_clean) <- crs
  
  #extract las from stems
  las_stems_rel <- stems[[6]]
  
  
  print("Calculating Tree Height...")
  #get height of tree
  i <- 1
  for (stemID in stems_sdf_clean$stemID){
    point <- stems_sdf_clean[i, ]
    circle <- st_buffer(point, dist = 1.5)
    ch <- max(terra::extract(chm, circle), na.rm = TRUE)
    stems_sdf_clean$ch[i] <- ch
    i <- i + 1
  }
  
  # get percentage of tree height
  stems_sdf_clean$h_dfl <- 0.6 * stems_sdf_clean$ch

  
  #calculate tree curvature
  value_df <- data.frame()
  i <- 1
  print("Calculating Stem Curvature...")
  
  for (ID_rel in stems_sdf_clean$stemID){
    #print(ID_rel)
    
    #get empty nn from above
    nn <- stems_sdf_clean$nn[stems_sdf_clean$stemID == ID_rel]
    h <- stems_sdf_clean$h_dfl[stems_sdf_clean$stemID == ID_rel]
    
    #filter stem single, transform data and normalise
    stem_single <- filter_poi(las_stems_rel, stemID == ID_rel & Z > 0.5 & Z < h)
    stem_df <- data.frame(cbind(stem_single$X, stem_single$Y, stem_single$Z))
    colnames(stem_df) <- c('X', 'Y', 'Z')
    stem_df$X_ref <- stem_df$X - min(stem_df$X)
    stem_df$Y_ref <- stem_df$Y - min(stem_df$Y)
    
    if (length(stem_df$X) > 10){
      #try with PCA aswell
      stem_df_pca <- stem_df[, c(3,4,5)]
      pca_res <- prcomp(stem_df_pca)
      growth_direction <- pca_res$rotation[, 1]
    }
    
    #aggreagte into 1cm slices using mean
    stem_df$Z_bin <- floor(stem_df$Z * 5) / 5
    stem_df_agg <- stem_df %>%
      group_by(Z_bin) %>%
      summarise(
        mean_Xr = median(X_ref, na.rm = TRUE),
        mean_Yr = median(Y_ref, na.rm = TRUE),
        count = n()
      )

    
    #######adapted version 2025
    #create spline fit for both dims
    spline_fit_X <- smooth.spline(stem_df_agg$Z_bin, stem_df_agg$mean_Xr, spar = 0.85)
    spline_fit_Y <- smooth.spline(stem_df_agg$Z_bin, stem_df_agg$mean_Yr, spar = 0.85)
    
    #len_splice <- length(spline_fit_Y$y)
    #print(len_splice)
    
    df_spline <- data.frame(cbind(spline_fit_X$x, spline_fit_X$y, spline_fit_Y$y))
    colnames(df_spline) <- c('Z_bin', 'X_spl', 'Y_spl')
    
    #calc mean for starting and end points points of line
    mean_first_Xr <- median(head(stem_df_agg$mean_Xr, 6))
    mean_last_Xr <- median(tail(stem_df_agg$mean_Xr, 6))
    mean_first_Yr <- median(head(stem_df_agg$mean_Yr, 6))
    mean_last_Yr <- median(tail(stem_df_agg$mean_Yr, 6))
    
    #create line fit to compare against
    lin_fit_df_x <- data.frame(y = c(mean_first_Xr, mean_last_Xr), x = c(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin)))
    line_fit <- lm(y ~ x, data = lin_fit_df_x)
    
    #get x and y from line_fit
    line_fit_x <- seq(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin), by = 0.1)
    lin_fitted_df_x <- data.frame(x = line_fit_x)
    lin_fitted_df_x$y <- predict(line_fit, newdata = lin_fitted_df_x)
    
    #create line fit to compare against
    lin_fit_df_y <- data.frame(y = c(mean_first_Yr, mean_last_Yr), x = c(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin)))
    line_fit <- lm(y ~ x, data = lin_fit_df_y)
    
    #get x and y from line_fit
    line_fit_y <- seq(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin), by = 0.1)
    lin_fitted_df_y <- data.frame(x = line_fit_y)
    lin_fitted_df_y$y <- predict(line_fit, newdata = lin_fitted_df_y)
    
    #create df of both lin fit and spline fit
    
    df_fits <- merge(stem_df_agg, lin_fitted_df_x, by.x = 'Z_bin', by.y = 'x')
    df_fits <- merge(df_fits, lin_fitted_df_y, by.x = 'Z_bin', by.y = 'x')
    
    colnames(df_fits) <- c('Z_bin', 'mean_Xr', 'meanYr', 'count', 'lin_x', 'lin_y')
    
    df_fits <- merge(df_fits, df_spline, by = 'Z_bin')

    
    #x DFL
    DFL_X <- sum(abs(df_fits$lin_x - df_fits$X_spl))/length(df_fits$lin_x)
    DFL_Y <- sum(abs(df_fits$lin_y - df_fits$Y_spl))/length(df_fits$lin_y)
    
    
    value_df <- rbind(value_df, c(ID_rel, DFL_X, DFL_Y, h / 0.6, st_coordinates(stems_sdf_clean[i,])[1], st_coordinates(stems_sdf_clean[i,])[2]))
    
    i <-  i + 1
  }
  
  colnames(value_df) <- c('stemID', 'DFL_X', 'DFL_Y', 'h', 'X', 'Y')
  value_df$DFL_sum <- value_df$DFL_X + value_df$DFL_Y
  
  #reorder cols
  #todo
  
  value_sdf <- st_as_sf(value_df, coords = c('X', 'Y'), crs = 32633)
  
  return(value_sdf)
}


##############################
#function for stem curve by DFL algo, required stem segmn, chm

pme_tree_curvature_plot_check <- function(stems, sID){
  require(dplyr)
  require(ggplot2)
  library(patchwork)
  #extract las from stems
  las_stems_rel <- filter_poi(stems[[6]], stemID == sID)
  
  
  #calculate tree curvature
  value_df <- data.frame()
  i <- 1
  #print("Calculating Stem Curvature...")
  
  for (ID_rel in stemID){
    #print(ID_rel
    
    #filter stem single, transform data and normalise
    stem_single <- filter_poi(las_stems_rel, stemID == stemID & Z > 0.5)
    stem_df <- data.frame(cbind(stem_single$X, stem_single$Y, stem_single$Z))
    colnames(stem_df) <- c('X', 'Y', 'Z')
    stem_df$X_ref <- stem_df$X - min(stem_df$X)
    stem_df$Y_ref <- stem_df$Y - min(stem_df$Y)
    
    if (length(stem_df$X) > 10){
      #try with PCA aswell
      stem_df_pca <- stem_df[, c(3,4,5)]
      pca_res <- prcomp(stem_df_pca)
      growth_direction <- pca_res$rotation[, 1]
    }
    
    #aggreagte into 1cm slices using mean
    stem_df$Z_bin <- floor(stem_df$Z * 5) / 5
    stem_df_agg <- stem_df %>%
      group_by(Z_bin) %>%
      summarise(
        mean_Xr = median(X_ref, na.rm = TRUE),
        mean_Yr = median(Y_ref, na.rm = TRUE),
        count = n()
      )
    
    
    #######adapted version 2025
    #create spline fit for both dims
    spline_fit_X <- smooth.spline(stem_df_agg$Z_bin, stem_df_agg$mean_Xr, spar = 0.85)
    spline_fit_Y <- smooth.spline(stem_df_agg$Z_bin, stem_df_agg$mean_Yr, spar = 0.85)
    
    #len_splice <- length(spline_fit_Y$y)
    #print(len_splice)
    
    df_spline <- data.frame(cbind(spline_fit_X$x, spline_fit_X$y, spline_fit_Y$y))
    colnames(df_spline) <- c('Z_bin', 'X_spl', 'Y_spl')
    
    #calc mean for starting and end points points of line
    mean_first_Xr <- median(head(stem_df_agg$mean_Xr, 6))
    mean_last_Xr <- median(tail(stem_df_agg$mean_Xr, 6))
    mean_first_Yr <- median(head(stem_df_agg$mean_Yr, 6))
    mean_last_Yr <- median(tail(stem_df_agg$mean_Yr, 6))
    
    #create line fit to compare against
    lin_fit_df_x <- data.frame(y = c(mean_first_Xr, mean_last_Xr), x = c(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin)))
    line_fit <- lm(y ~ x, data = lin_fit_df_x)
    
    #get x and y from line_fit
    line_fit_x <- seq(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin), by = 0.1)
    lin_fitted_df_x <- data.frame(x = line_fit_x)
    lin_fitted_df_x$y <- predict(line_fit, newdata = lin_fitted_df_x)
    
    #create line fit to compare against
    lin_fit_df_y <- data.frame(y = c(mean_first_Yr, mean_last_Yr), x = c(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin)))
    line_fit <- lm(y ~ x, data = lin_fit_df_y)
    
    #get x and y from line_fit
    line_fit_y <- seq(min(stem_df_agg$Z_bin), max(stem_df_agg$Z_bin), by = 0.1)
    lin_fitted_df_y <- data.frame(x = line_fit_y)
    lin_fitted_df_y$y <- predict(line_fit, newdata = lin_fitted_df_y)
    
    #create df of both lin fit and spline fit
    
    df_fits <- merge(stem_df_agg, lin_fitted_df_x, by.x = 'Z_bin', by.y = 'x')
    df_fits <- merge(df_fits, lin_fitted_df_y, by.x = 'Z_bin', by.y = 'x')
    
    colnames(df_fits) <- c('Z_bin', 'mean_Xr', 'meanYr', 'count', 'lin_x', 'lin_y')
    
    df_fits <- merge(df_fits, df_spline, by = 'Z_bin')
    
    
    
    df_fits <- df_fits[order(df_fits$Z_bin), ]
    #x DFL
    DFL_X <- sum(abs(df_fits$lin_x - df_fits$X_spl))/length(df_fits$lin_x)
    DFL_Y <- sum(abs(df_fits$lin_y - df_fits$Y_spl))/length(df_fits$lin_y)
    
    plot_title <- paste("DFL_X:", round(DFL_X, 4), "| DFL_Y:", round(DFL_Y, 4), round((DFL_Y + DFL_X), 4))
    
    # Build the plot
    p_xz <- ggplot() +
      geom_point(data = stem_df, aes(x = X_ref, y = Z), color = "black", alpha = 0.5) +
      geom_line(data = df_fits, aes(x = lin_x, y = Z_bin), color = "blue", size = 1, linetype = "dashed") +
      geom_point(data = df_fits, aes(x = X_spl, y = Z_bin), color = "red", size = 3) +
      labs(title = "X-Z Plot", x = "X", y = "Z") + coord_fixed() + 
      theme_minimal()
    
    p_yz <- ggplot() +
      geom_point(data = stem_df, aes(x = Y_ref, y = Z), color = "black", alpha = 0.5) +
      geom_line(data = df_fits, aes(x = lin_y, y = Z_bin), color = "blue", size = 1, linetype = "dashed") +
      geom_point(data = df_fits, aes(x = Y_spl, y = Z_bin), color = "red", size = 3) +
      labs(title = "Y-Z Plot", x = "Y", y = "Z") + coord_fixed() + 
      theme_minimal()
    
    
    combined_plot <- p_xz + p_yz + plot_layout(ncol = 2)  +
      plot_annotation(title = plot_title)
    
    
    i <-  i + 1
  }
  return(combined_plot)
}
