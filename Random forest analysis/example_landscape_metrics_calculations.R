# Bring in contact rates
contact.est.slow.decay <- readRDS("Example simulation/Inputs/contact_rates_annual_lone_oaks.RDS")
contact.est.fast.decay <- readRDS("Example simulation/Inputs/contact_rates_annual_CP.RDS")

slow_kernel_df <-  tibble(distance_km = seq(0.001, 5, by = 0.001), 
                          contact = predict(contact.est.slow.decay, 
                                            newdata = tibble(distance = seq(0.001, 5, by = 0.001)), 
                                            re.form = NA, type = "response")) %>% 
  filter(contact >= 0.001)

fast_kernel_df <- tibble(distance_km = seq(0.001, 3, by = 0.001), 
                         contact = predict(contact.est.fast.decay, 
                                           newdata = tibble(distance = seq(0.001, 3, by = 0.001)), 
                                           re.form = NA, type = "response")) %>% 
  filter(contact >= 0.001)


extract_multiscale_metrics <- function(
    file,
    grid,
    nlcd_file_path,
    contact_df,
    pixel_distance_df,
    study_area,
    season,
    contact_category = NA
){
  
  library(dplyr)
  library(sf)
  library(tibble)
  library(terra)
  library(matrixStats)
  library(landscapemetrics)
  library(purrr)
  library(tidyterra)
  
  input <- readRDS(file)
  
  origin_df <- input$origin_probs
  selection_df <- input$selection_probs
  
  #--------------------------------------------------
  # NLCD
  #--------------------------------------------------
  
  nlcd_rast <- terra::rast(nlcd_file_path)
  
  rcl_mat <- matrix(c(
    11, 1, 21, 2, 22, 2, 23, 3, 24, 3, 31, 3,
    41, 4, 42, 4, 43, 4, 52, 4,
    71, 5, 81, 5, 82, 5,
    90, 6, 95, 6
  ), ncol = 2, byrow = TRUE)
  
  reclass_rast <- terra::classify(nlcd_rast, rcl_mat)
  
  names(reclass_rast) <- "landcover"
  
  reclass_rast <- terra::crop(reclass_rast,grid,mask = TRUE) %>%
    tidyterra::mutate(landcover = factor(landcover))
  
  #reclass_rast <- terra::crop(reclass_rast, grid, mask = TRUE)
  
  safe_metric <- function(expr){
    tryCatch(expr, error = function(e) NA_real_)
  }
  
  #--------------------------------------------------
  # CELL SCALE
  #--------------------------------------------------
  
  origin_mat <- as.matrix(origin_df[, -1])
  selection_mat <- as.matrix(selection_df[, -1])
  
  cell_metrics <- tibble(
    grid_id = origin_df$grid_id,
    
    mean_origin_cell = rowMeans(origin_mat, na.rm = TRUE),
    sd_deer_origin_cell = apply(origin_mat, 1, sd, na.rm = TRUE),
    
    mean_selection_cell = rowMeans(selection_mat, na.rm = TRUE),
    sd_deer_selection_cell = apply(selection_mat, 1, sd, na.rm = TRUE)
  )
  
  #--------------------------------------------------
  # LANDSCAPE SCALE
  #--------------------------------------------------
  
  landscape_metrics <- tibble(
    mean_origin_landscape =
      mean(cell_metrics$mean_origin_cell, na.rm = TRUE),
    
    sd_origin_landscape =
      sd(cell_metrics$mean_origin_cell, na.rm = TRUE),
    
    sd_deer_origin_landscape =
      median(cell_metrics$sd_deer_origin_cell, na.rm = TRUE),
    
    mean_selection_landscape =
      mean(cell_metrics$mean_selection_cell, na.rm = TRUE),
    
    sd_selection_landscape =
      sd(cell_metrics$mean_selection_cell, na.rm = TRUE),
    
    sd_deer_selection_landscape =
      median(cell_metrics$sd_deer_selection_cell, na.rm = TRUE),
    
    agg_index_landscape =
      safe_metric(lsm_l_ai(reclass_rast, directions = 8)$value),
    
    cohesion_landscape =
      safe_metric(lsm_l_cohesion(reclass_rast, directions = 8)$value),
    
    division_landscape =
      safe_metric(lsm_l_division(reclass_rast, directions = 8)$value)
  )
  
  #--------------------------------------------------
  # CELL FRAGMENTATION
  #--------------------------------------------------
  
  cell_frag <- purrr::map_dfr(grid$grid_id, function(id){
    
    g <- grid |> tidyterra::filter(grid_id == id)
    
    r <- tryCatch(
      terra::crop(reclass_rast, g, mask = TRUE),
      error = function(e) NULL
    )
    
    if(is.null(r)){
      return(tibble(
        grid_id = id,
        agg_index_cell = NA_real_,
        cohesion_cell = NA_real_,
        division_cell = NA_real_
      ))
    }
    
    tibble(
      grid_id = id,
      agg_index_cell =
        safe_metric(lsm_l_ai(r, directions = 8)$value),
      cohesion_cell =
        safe_metric(lsm_l_cohesion(r, directions = 8)$value),
      division_cell =
        safe_metric(lsm_l_division(r, directions = 8)$value)
    )
  })
  
  #--------------------------------------------------
  # KERNEL LOOKUP
  #--------------------------------------------------
  
  max_dist <- max(contact_df$distance_km)
  
  cell_lookup <- pixel_distance_df %>%
    dplyr::select(pixel_1, pixel_2, distance) %>%
    mutate(distance_km = round(distance / 1000, 3)) %>%
    filter(distance_km <= max_dist) %>%
    dplyr::select(pixel_1, pixel_2) %>%
    distinct()
  
  #--------------------------------------------------
  # KERNEL METRICS
  #--------------------------------------------------
  
  kernel_metrics <- purrr::map_dfr(grid$grid_id, function(cell_id){
    
    grid_ids <- cell_lookup %>%
      filter(pixel_1 == cell_id) %>%
      pull(pixel_2)
    
    sel <- selection_df %>%
      filter(grid_id %in% grid_ids) %>%
      dplyr::select(-1) %>%
      as.matrix()
    
    sel <- sel[rowSums(sel != 0) > 0, , drop = FALSE]
    
    if(nrow(sel) > 0){
      
      mean_selection_kernel <- mean(rowMeans(sel), na.rm = TRUE)
      sd_selection_kernel <- sd(rowMeans(sel), na.rm = TRUE)
      sd_deer_selection_kernel <- median(
        apply(sel, 1, sd, na.rm = TRUE),
        na.rm = TRUE
      )
      
    } else {
      
      mean_selection_kernel <- NA_real_
      sd_selection_kernel <- NA_real_
      sd_deer_selection_kernel <- NA_real_
    }
    
    org <- origin_df %>%
      filter(grid_id %in% grid_ids) %>%
      dplyr::select(-1) %>%
      as.matrix()
    
    org <- org[rowSums(org != 0) > 0, , drop = FALSE]
    
    if(nrow(org) > 0){
      
      mean_origin_kernel <- mean(rowMeans(org), na.rm = TRUE)
      sd_origin_kernel <- sd(rowMeans(org), na.rm = TRUE)
      sd_deer_origin_kernel <- median(
        apply(org, 1, sd, na.rm = TRUE),
        na.rm = TRUE
      )
      
    } else {
      
      mean_origin_kernel <- NA_real_
      sd_origin_kernel <- NA_real_
      sd_deer_origin_kernel <- NA_real_
    }
    
    gk <- grid |> tidyterra::filter(grid_id %in% grid_ids)
    
    rk <- tryCatch(
      terra::crop(reclass_rast, gk, mask = TRUE),
      error = function(e) NULL
    )
    
    tibble(
      grid_id = cell_id,
      
      mean_origin_kernel = mean_origin_kernel,
      sd_origin_kernel = sd_origin_kernel,
      sd_deer_origin_kernel = sd_deer_origin_kernel,
      
      mean_selection_kernel = mean_selection_kernel,
      sd_selection_kernel = sd_selection_kernel,
      sd_deer_selection_kernel = sd_deer_selection_kernel,
      
      agg_index_kernel = if(is.null(rk)) NA_real_
      else safe_metric(lsm_l_ai(rk, directions = 8)$value),
      
      cohesion_kernel = if(is.null(rk)) NA_real_
      else safe_metric(lsm_l_cohesion(rk, directions = 8)$value),
      
      division_kernel = if(is.null(rk)) NA_real_
      else safe_metric(lsm_l_division(rk, directions = 8)$value)
    )
  })
  
  #--------------------------------------------------
  # FINAL OUTPUT
  #--------------------------------------------------
  
  output <- cell_metrics %>%
    left_join(cell_frag, by = "grid_id") %>%
    left_join(kernel_metrics, by = "grid_id") %>%
    mutate(
      study_area = study_area,
      season = season,
      contact_category = contact_category
    ) %>%
    dplyr::bind_cols(landscape_metrics)
  
  output
}

##ECRP ####

ECRP_fall_slow <- extract_multiscale_metrics(file = "Example simulation/Inputs/ECRP_sim_input_fall.RDS", 
                                             study_area = "ECRP", 
                                             season = "fall", 
                                             grid = readRDS("Example simulation/Inputs/ECRP_grid.RDS"),
                                             nlcd_file_path = "Example simulation/Inputs/ECRP_NLCD.RDS",
                                             contact_df = slow_kernel_df,
                                             pixel_distance_df = readRDS("Example simulation/Inputs/ECRP_grid_dist.RDS"),
                                             contact_category = "slow")

ECRP_fall_fast <- extract_multiscale_metrics(file = "~/ER_Data_Files/simulation_inputs/ECRP/ECRP_sim_input_fall.RDS", 
                                             study_area = "ECRP", 
                                             season = "fall", 
                                             grid = readRDS("~/ER_Data_Files/study_areas/ECRP/ECRP_grid.RDS"),
                                             nlcd_file_path = "~/ER_Data_Files/land_cover/ECRP/nlcd_2021_land_cover_metro_area_clipped.tif",
                                             contact_df = fast_kernel_df,
                                             pixel_distance_df = readRDS("~/ER_Data_Files/grid_cell_distances/ECRP_grid_dist.RDS"),
                                             contact_category = "fast")

