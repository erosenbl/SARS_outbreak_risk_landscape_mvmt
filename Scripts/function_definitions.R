# Script to load all packages and functions for outbreak analysis

#Load all packages
library(tidyverse) 
library(lubridate) 
library(sf) 
library(raster) 
library(amt) 
library(data.table) 
library(wildlifeDI) 
library(move2) 
library(glmmTMB) 
library(terra) 
library(lme4) 
library(extraDistr) 
library(whitetailedSIRS) 
library(adehabitatHR) 
library(tidyterra) 
library(drc) 
library(parallel) 
library(forcats) 
library(Hmisc) 
library(units) 
library(IndRSA) 
library(DescTools) 
library(future.apply) 
library(landscapemetrics) 
library(matrixStats) 

####
# GPS_fix_prep: Helper function to prepare collar locations from either csv or 
# shp, as study area datasets have slightly different formatting
GPS_fix_prep <- function(directory = NULL, vector_object = NULL, start, end){
  require(tidyverse)
  require(lubridate)
  
  if(!is.null(directory)){
    # Get the list of all CSV files in the directory
    collar_data <- readRDS(file = directory)}
  
  if(!is.null(vector_object)){
    collar_data <- vector_object
  }
  
  #Filter by season
  # Define date range for filtering - FALL
  start_date <- as.Date(start)
  end_date <- as.Date(end)
  
  #Filter by start and end dates
  collar_data_filtered <- collar_data %>% 
    filter(date_time >= start_date & date_time <= end_date)
  
  collar_data_filtered <- as.data.frame(collar_data_filtered)
  
  return(collar_data_filtered)
}

####
# collar.meta.prep: Helper function to prepare collar metadata from csv file, 
# as study area datasets have slightly different formatting
collar.meta.prep <- function(filepath, Site_name, capture_field, capture_status, start_field, end_field){
  require(tidyverse)
  #Bring in metadata to assign sex and age classes
  collar.meta <- read.csv(filepath, na.strings = c("", "NA", "<NA>"))
  
  collar.meta <- collar.meta %>% 
    filter(.data[[capture_field]] == capture_status)
  
  collar.meta <- collar.meta %>% 
    filter(!is.na(animal_id), site_name == Site_name) %>% 
    dplyr::select(., animal_id,sex, age, start_field, end_field) %>% 
    mutate(age = case_when(age %in% c("fawn", "YoungOfYear", "Young of year")~ "Fawn",
                           age == "yearling" ~ "Yearling",
                           age == "adult" ~ "Adult",
                           TRUE ~ age)) %>% 
    rename(., "Sex" = "sex", "Age_Class" = "age", 
           "start_deployment" = start_field, 
           "end_deployment" = end_field) %>% 
    mutate(., start_deployment = as.Date(start_deployment, format = "%m/%d/%Y"),
           end_deployment = as.Date(end_deployment, format = "%m/%d/%Y")) %>% 
    distinct()
  
  return(collar.meta)
}

####
# collar.meta.prep: Prepare collar metadata from csv file, if metafile data is 
# nested with collar location data
meta_from_collars <- function(directory, site_name, columns_to_keep){
  # Get the list of all CSV files in the directory
  file_list <- list.files(path = directory, pattern = "*.csv", 
                          full.names = TRUE)
  
  # Function to read each CSV file
  read_csv_files <- function(file) {
    read.csv(file, stringsAsFactors = FALSE)
  }
  
  # Read all CSV files into a list of data frames
  data_list <- lapply(file_list, read_csv_files)
  
  # Renaming the list based on the value in the first row, second column of 
  # each data frame
  names(data_list) <- sapply(data_list, function(df) df[1, 2])
  
  collar.meta <- rbindlist(data_list) %>% 
    dplyr::select(all_of(columns_to_keep)) %>% 
    filter(site == site_name) %>% 
    distinct() %>% 
    mutate(sex = case_when(sex %in% c("f","F","female","Female") ~ "Female",
                           sex %in% c("m", "M", "male", "Male") ~ "Male"))
  return(collar.meta)
}


####
# centroid_dist: Calculate pairwise distances between collared animals' 
# centroids

centroid_dist <- function(collar_locations, centroid_crs){
  
  require(tidyverse)
  require(lubridate)
  require(sf)
  require(raster)
  
  loc_by_id <- collar_locations %>% 
    dplyr::select(animal_id = deployid, date_time, x, y) %>% 
    merge(., collar.meta, by = c("deploy_id" = "animal_id")) %>% 
    st_as_sf(., coords = 3:4, crs = centroid_crs)
  
  
  ##STEP 1: Identify centroid grid cell for each animal
  # Generate hr centroid for each animal, and add the id column to get unique 
  # home ranges [1]
  CR_deer_ud <- adehabitatHR::kernelUD(as(loc_by_id, "Spatial")[1], grid = 500, 
                                       extent = 0.5)
  
  #Extract home ranges for each animal
  hr_vol <- adehabitatHR::getvolumeUD(CR_deer_ud, standardize = T)
  
  ud_rasters <- lapply(1:length(hr_vol), function(i) {
    raster(as(hr_vol[[i]], "SpatialPixelsDataFrame"))
  })
  
  # Assign animal IDs as names
  names(ud_rasters) <- names(hr_vol)  
  
  # Get coordinates of all pixels and extract UD values
  centroids <- lapply(ud_rasters, function(r) {
    coords <- xyFromCell(r, 1:ncell(r))  
    values <- getValues(r)
    
    centroid <- coords[which(values==min(values)),]
    return(centroid)
  })
  
  # Convert to a data frame for easy reference
  centroid_df <- as_tibble(do.call(rbind, centroids)) %>% 
    mutate(animal_id = names(ud_rasters))
  
  # Step 2: Convert centroids to an sf object for spatial distance calculation, 
  #with the option to assigne CRS if not lat/long
  centroids_sf <- st_as_sf(centroid_df, coords = c("x", "y"), crs = 
                             centroid_crs)
  
  # Step 3: Calculate pairwise distances between centroids
  dist_matrix <- st_distance(centroids_sf)
  
  # Step 4: Convert distance matrix to a tidy data frame
  dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
  names(dist_df) <- c("deployid1", "deployid2", "distance")
  
  # Step 5: Merge with centroids_df to include deployid information
  dist_df <- dist_df %>%
    mutate(
      deployid1 = centroid_df$animal_id[as.integer(deployid1)],
      deployid2 = centroid_df$animal_id[as.integer(deployid2)],
      distance = as.numeric(distance)
    ) %>%
    filter(deployid1 != deployid2) %>% 
    rowwise() %>%
    mutate(
      deployid_sorted1 = min(deployid1, deployid2),
      deployid_sorted2 = max(deployid1, deployid2)
    ) %>%
    ungroup() %>%
    dplyr::select(deployid = deployid_sorted1, contact_id = deployid_sorted2, 
                  distance) %>%
    distinct(deployid, contact_id, .keep_all = TRUE)  #
  
  # View the result
  return(dist_df)
}

####
# daily_centroid_dist: Calculate pairwise distances between collared animals' 
# centroids on a daily basis

daily_centroid_dist <- function(animal_data, collar_meta, centroid_crs,
                                grid_size = 30, extent_margin = 500,
                                min_points_per_day = 5,
                                num_cores = 2) {
  require(tidyverse)
  require(lubridate)
  require(sf)
  require(raster)
  require(adehabitatHR)
  require(parallel)
  
  # Step 1: Prepare location data
  loc_by_id <- animal_data %>%
    left_join(collar_meta, by = c("deployid" = "animal_id")) %>%
    mutate(date = as.Date(date_time)) %>%
    rename("animal_id" = "deployid") %>%
    st_as_sf(coords = c("x", "y"), crs = centroid_crs)
  
  # Step 2: Create fixed raster grid
  bbox_all <- st_bbox(loc_by_id)
  raster_template <- raster(
    xmn = bbox_all$xmin - extent_margin,
    xmx = bbox_all$xmax + extent_margin,
    ymn = bbox_all$ymin - extent_margin,
    ymx = bbox_all$ymax + extent_margin,
    resolution = grid_size
  )
  grid_sp <- as(raster_template, "SpatialPixels")
  
  # Step 3: Generate full date sequence
  
  date_range <- seq(min(loc_by_id$date), max(loc_by_id$date), by = "day")
  animals <- unique(loc_by_id$animal_id)
  full_dates <- expand.grid(animal_id = animals, date = date_range)
  
  # Step 4: Define helper function for daily centroid
  calculate_daily_centroid <- function(day_data, animal_id, date) {
    if (nrow(day_data) < min_points_per_day) {
      return(tibble(animal_id = animal_id, date = date, x = NA_real_, 
                    y = NA_real_))
    }
    
    # convert to Spatial for adehabitatHR
    day_data$animal_id <- animal_id
    sp_obj <- as(day_data[, c("animal_id","geometry")], "Spatial")
    ud <- kernelUD(sp_obj, grid = grid_sp)
    vol <- getvolumeUD(ud[[1]], standardize = TRUE)
    
    r <- raster(as(vol, "SpatialPixelsDataFrame"))
    coords <- xyFromCell(r, which.min(getValues(r)))
    
    tibble(animal_id = animal_id, date = date, x = coords[1], y = coords[2])
  }
  
  loc_split_tbl <- loc_by_id %>%
    group_by(animal_id, date) %>%
    nest()
  
  centroid_results <- mclapply(seq_len(nrow(loc_split_tbl)), function(i) {
    df <- loc_split_tbl$data[[i]]
    key_animal <- loc_split_tbl$animal_id[i]
    key_date <- loc_split_tbl$date[i]
    
    tryCatch(
      calculate_daily_centroid(df, key_animal, key_date),
      error = function(e) {
        message("Centroid error for animal: ", key_animal, " date: ", key_date, 
                " -> ", e$message)
        return(tibble(animal_id = key_animal, date = key_date, x = NA_real_, 
                      y = NA_real_))
      }
    )
  }, mc.cores = num_cores)
  
  daily_centroids <- bind_rows(centroid_results)
  
  # Ensure full animal-date coverage
  daily_centroids_full <- full_dates %>%
    left_join(daily_centroids, by = c("animal_id", "date"))
  
  # Step 6: Split by date for pairwise distance calculation
  
  date_list <- split(daily_centroids_full, daily_centroids_full$date)
  
  process_day <- function(day_data) {
    current_date <- unique(day_data$date)
    
    # filter out days with no centroid positions
    day_data <- day_data %>% filter(!is.na(x) & !is.na(y))
    if (nrow(day_data) < 2) {
      return(tibble(
        date = current_date,
        deployid = character(),
        contact_id = character(),
        distance = numeric()
      ))
    }
    
    day_sf <- st_as_sf(day_data, coords = c("x", "y"), crs = centroid_crs)
    dist_matrix <- st_distance(day_sf)
    
    dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
    names(dist_df) <- c("deployid1", "deployid2", "distance")
    
    dist_df <- dist_df %>%
      mutate(
        deployid1 = day_data$animal_id[as.integer(deployid1)],
        deployid2 = day_data$animal_id[as.integer(deployid2)],
        distance = as.numeric(distance)
      ) %>%
      filter(deployid1 != deployid2) %>%
      rowwise() %>%
      mutate(
        deployid_sorted1 = min(deployid1, deployid2),
        deployid_sorted2 = max(deployid1, deployid2)
      ) %>%
      ungroup() %>%
      dplyr::select(
        deployid = deployid_sorted1,
        contact_id = deployid_sorted2,
        distance
      ) %>%
      distinct(deployid, contact_id, .keep_all = TRUE) %>%
      mutate(date = current_date) %>%
      dplyr::select(date, deployid, contact_id, distance)
    
    return(dist_df)
  }
  
  # Step 7: Run distance calculation in parallel
  results_list <- mclapply(date_list, function(d) {
    tryCatch(
      process_day(d),
      error = function(e) {
        message("Distance error on date: ", unique(d$date), " -> ", e$message)
        return(NULL)
      }
    )
  }, mc.cores = num_cores)
  
  # Combine results
  pairwise_daily_distances <- bind_rows(results_list)
  
  return(pairwise_daily_distances)
}

####
# daily_step_length: Summarize daily step lengths for each radio-collared animal
# in a study area

daily_step_length <- function(data, CRS = "epsg:32615"){
  require(tidyverse)
  require(amt)
  
  #Nest data, create track, and calculate various movement metrics
  mvmt_metrics <- data  %>%
    nest(data = -deployid) %>%
    mutate(track = map(data, ~ make_track(tbl = .x, .x = x, .y = y,
                                          .t = date_time, crs = CRS))) %>%
    mutate(dat.resample = map(track, ~track_resample(., rate = days(1), 
                                                     tolerance = hours(1)))) %>%
    mutate(dir_abs = map(dat.resample, direction_abs, full_circle = T, 
                         zero = "N", clockwise = T),
           dir_rel = map(dat.resample, direction_rel),
           sl = map(dat.resample, step_lengths)) %>%
    unnest(cols = c(dir_abs, dir_rel, sl)) %>%
    dplyr::select(deployid, sl, dir_abs, dir_rel)
  
  return(mvmt_metrics)
}

####
# extract_movement: Fits gamma regression for daily step length as a function of
# sex and age class combinations, either calculating median daily displacement 
# (result = "median.mvmt") or gamma shape parameters 
# (result = "mvmt.dist.shape").

extract_movement <- function(daily_steps, collar_metadata, result){
  #Load required libraries
  require(tidyverse)
  require(fitdistrplus)
  require(lme4)
  
  if(!(result %in% c("median.mvmt","mvmt.dist.shape"))) {
    print("Please set result argument to either 'median.mvmt' or 
          'mvmt.dist.shape'.")
    stop()
  }
  
  #Joining daily movement data with animal sex and age class
  mvmt <- daily_steps %>% 
    merge(collar_metadata, by.x = "deployid", by.y = "animal_id") %>% 
    filter(!is.na(sl), sl>0) %>% 
    mutate(Sex = case_when(Sex == "F" ~ "Female",
                           Sex == "M" ~ "Male",
                           TRUE ~ Sex)) %>% 
    mutate(Sex = factor(Sex),
           Age_Class = factor(Age_Class)) %>% 
    filter(Age_Class %in% c("Adult"))#,"Yearling"))
  
  #Fit gamma regression with random effect of individual
  mvmt.mod <- glmer(sl~Sex+(1|deployid), data = mvmt, 
                    family = Gamma(link = "log"))
  
  new.data <- expand.grid(Sex = c("Female", "Male"), Age_Class = c("Adult"))
  
  pred <- predict(mvmt.mod, newdata = new.data, se.fit = T, re.form = ~ 0)
  pred.df <- cbind(new.data,est = exp(pred$fit), 
                   LCL95 = exp(pred$fit -1.96*pred$se.fit), 
                   UCL95 = exp(pred$fit + 1.96*pred$se.fit))
  
  if(result == "median.mvmt") {return(pred.df)}
  
  if(result == "mvmt.dist.shape")
  {mvmt.ad.yr <- mvmt %>%
    filter(., Age_Class %in% c("Adult")) %>% 
    mutate(Age_Sex = case_when(Age_Class == "Adult" & 
                                 Sex == "Female"~"Adult female",
      Age_Class == "Adult" & Sex == "Male"~"Adult male"))
  
  AF <- mvmt.ad.yr %>% filter(., Age_Sex == "Adult female")
  AF.fit.gamma <- fitdist(AF$sl, "gamma")
  
  AM <- mvmt.ad.yr %>% filter(., Age_Sex == "Adult male")
  AM.fit.gamma <- fitdist(AM$sl, "gamma", method = "mme")
  
  #Combine shape parameter estimates
  mvmt.shape <- tibble(Age_Sex = c("Adult female","Adult male"), 
                       shape = c(AF.fit.gamma$estimate[1], 
                                 AM.fit.gamma$estimate[1]),
                       rate = c(AF.fit.gamma$estimate[2], 
                                AM.fit.gamma$estimate[2]))
  
  return(mvmt.shape)}
  
}

####
# identify_contacts: A function to identify proximity events
# Modified from wildlifeDI package 
# (https://github.com/jedalong/wildlifeDI/blob/master/R/conProcess.R)

library(parallel)
library(dplyr)
library(units)

identify_contacts <- function(traj, traj2, dc = 0, tc = 0, GetSim = TRUE, fixid,
                              num_cores = 2) {
  
  # global vars
  id1 <- NULL
  dist <- NULL
  fixid1 <- NULL
  
  # Unit control
  units(tc) <- as_units("s")
  
  # Combine trajectories and identify overlap pairs
  if (missing(traj2)) {
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO == TRUE, ]
    mtraj <- traj
  } else {
    pairs <- checkTO(traj, traj2)
    pairs <- pairs[pairs$TO == TRUE, ]
    if (st_crs(traj2) != st_crs(traj)) {
      traj2 <- st_transform(traj2, crs = st_crs(traj))
    }
    mtraj <- mt_stack(traj, traj2, track_combine = 'check_unique')
  }
  
  # Set up rownames
  if (missing(fixid)) {
    id_mtraj <- as.character(mt_track_id(mtraj))
    row.names(mtraj) <- paste0(id_mtraj, '_', stats::ave(id_mtraj, id_mtraj, 
                                                         FUN = seq_along))
    id_traj <- as.character(mt_track_id(traj))
    row.names(traj) <- paste0(id_traj, '_', stats::ave(id_traj, id_traj, 
                                                       FUN = seq_along))
  } else {
    row.names(mtraj) <- mtraj[[fixid]]
  }
  
  n.pairs <- nrow(pairs)
  
  # Parallel loop with safe error handling
  results <- mclapply(1:n.pairs, function(i) {
    tryCatch({
      # Subset tracks
      traja <- mtraj[mt_track_id(mtraj) == pairs$ID1[i], ]
      trajb <- mtraj[mt_track_id(mtraj) == pairs$ID2[i], ]
      
      # Skip if either track is empty
      if (nrow(traja) == 0 || nrow(trajb) == 0) {
        message("Skipping empty track at pair ", i, ": ", pairs$ID1[i], " vs ", 
                pairs$ID2[i])
        return(NULL)
      }
      
      if (GetSim) {
        trajs <- GetSimultaneous(traja, trajb, tc)
        tr1 <- trajs[mt_track_id(trajs) == pairs$ID1[i], ]
        tr2 <- trajs[mt_track_id(trajs) == pairs$ID2[i], ]
        
        # Skip if GetSimultaneous produced mismatched data
        if (nrow(tr1) == 0 || nrow(tr2) == 0) {
          message("No simultaneous fixes within tc for pair ", i, ": ", 
                  pairs$ID1[i], " vs ", pairs$ID2[i])
          return(NULL)
        }
        
        proxdf <- data.frame(
          id1 = pairs$ID1[i],
          id2 = pairs$ID2[i],
          fixid1 = row.names(tr1),
          fixid2 = row.names(tr2),
          t1 = mt_time(tr1),
          t2 = mt_time(tr2),
          dist = st_distance(tr1, tr2, by_element = TRUE),
          difftime = as.numeric(abs(mt_time(tr1) - mt_time(tr2)))
        )
        units(proxdf$difftime) <- as_units('s')
        units(dc) <- units(proxdf$dist)
        
        # Filter based on distance and time thresholds
        proxdf <- proxdf[proxdf$dist < dc & proxdf$difftime < tc, ]
        
        # If no rows left, return NULL
        if (nrow(proxdf) == 0) return(NULL)
        
      } else {
        # Non-GetSim approach
        dM <- st_distance(traja, trajb)
        tM <- abs(-outer(as.numeric(mt_time(traja)), as.numeric(mt_time(trajb)),
                         '-'))
        
        units(tM) <- as_units('s')
        units(dc) <- units(dM)
        
        rownames(tM) <- rownames(traja)
        colnames(tM) <- rownames(trajb)
        
        ind <- which(dM < dc & tM < tc, arr.ind = TRUE)
        
        if (length(ind) == 0) return(NULL)
        
        rnm1 <- rownames(tM)[ind[, 1]]
        rnm2 <- colnames(tM)[ind[, 2]]
        
        proxdf <- data.frame(
          id1 = pairs$ID1[i],
          id2 = pairs$ID2[i],
          fixid1 = rnm1,
          fixid2 = rnm2,
          t1 = mt_time(traja)[ind[, 1]],
          t2 = mt_time(trajb)[ind[, 2]],
          dist = dM[ind],
          difftime = tM[ind]
        )
      }
      
      return(proxdf)
      
    }, error = function(e) {
      # Catch errors and safely return NULL
      message("Error at pair ", i, ": ", pairs$ID1[i], " vs ", pairs$ID2[i], 
              " -> ", e$message)
      return(NULL)
    })
    
  }, mc.cores = num_cores)
  condf <- bind_rows(results)
  return(condf)
}

####
# calc_contact_rates_continuous_gamma: A function to calculate the daily 
# probability of contact and the frequency of contacts in a given day when 
# contact occurs.

calc_contact_rates_continuous_gamma <- function(collar_locations, 
                                contact_crs = '+proj=utm +zone=15 ellps=WGS84', 
                                fix_animal_field = "deployid",
                                time_field = "date_time", 
                                distance_tol = 25, 
                                time_tol = 15*60, 
                                collar_metadata, 
                                collar_metadata_animal_field = "animal_id", 
                                centroid_distances, 
                                distance_interval = "daily", 
                                median_mvmt, 
                                start.date, 
                                end.date, 
                                num_cores = 1) {
  
  require(tidyverse)
  require(lubridate)
  require(sf)
  require(data.table)
  require(wildlifeDI)
  require(move2)
  require(glmmTMB)
  
  ## Convert locations to sf
  print("Setting up locations as sf object")
  sf_locs <- sf::st_as_sf(collar_locations, coords = c("x", "y")) %>% 
    sf::st_set_crs(contact_crs)
  
  move_df <- move2::mt_as_move2(sf_locs, time_column = time_field, 
                                track_id_column = fix_animal_field)
  
  ## Generate contacts using wildlifeDI
  print("Generating contact events. This will take 10-20 minutes...")
  move_df_cons <- identify_contacts(traj = move_df, dc = distance_tol, 
                                    tc = time_tol, num_cores = num_cores)
  
  # STEP 1: Create full set of possible dyads and dates
  ### NEW: Get all unique animals
  all_animals <- unique(collar_locations[[fix_animal_field]])
  
  ###Create all possible unique dyads
  all_dyads <- t(combn(all_animals, 2)) %>%
    as_tibble() %>%
    rename(deployid = V1, contact_id = V2)
  
  ###Create full date range
  date_range <- seq(as.Date(start.date), as.Date(end.date), by = "day")
  
  ###Expand grid to include every dyad × every date
  dyad_date_grid <- expand_grid(all_dyads, date = date_range)
  
  ###Add centroid distances for each dyad
  dyad_date_grid <- dyad_date_grid %>%
    left_join(centroid_distances, by = c("deployid", "contact_id", "date"))
  
  # STEP 2: Summarize identified contacts by day
  daily_contact_summary <- move_df_cons %>%
    mutate(
      date = as.Date(t1),
      # sort IDs to avoid duplicate pairings
      deployid_sorted = pmin(id1, id2),
      contact_id_sorted = pmax(id1, id2)
    ) %>%
    group_by(deployid_sorted, contact_id_sorted, date) %>%
    summarise(
      daily_contact = as.integer(n() >= 1),  # at least one contact that day
      contact_count = n(),                   # total count that day
      .groups = "drop"
    ) %>%
    rename(deployid = deployid_sorted, contact_id = contact_id_sorted)

  if(distance_interval == "season") {
    
    ###Merge full grid with observed contacts
    complete_data <- dyad_date_grid %>%
      left_join(daily_contact_summary, 
                by = c("deployid", "contact_id", "date")) %>%
      mutate(
        daily_contact = replace_na(daily_contact, 0),
        contact_count = replace_na(contact_count, 0)
      )
    
    # Filter to only dates within deployment period
    complete_contact_summary_Age_Sex <- complete_data %>% 
      merge(., collar_metadata, by.x = "deployid", by.y = "animal_id") %>% 
      dplyr::select(-Sex, -Age_Class) %>%
      dplyr::rename(start_deployment_deploy = start_deployment, 
                    end_deployment_deploy = end_deployment) %>%
      merge(., collar_metadata, by.x = "contact_id", by.y = "animal_id") %>% 
      dplyr::select(-Sex, -Age_Class) %>%
      dplyr::rename(start_deployment_contact = start_deployment, 
                    end_deployment_contact = end_deployment) %>%
      filter(date >= start_deployment_deploy &
               (is.na(end_deployment_deploy) | date <= end_deployment_deploy),
             date >= start_deployment_contact &
               (is.na(end_deployment_contact) | 
                  date <= end_deployment_contact)) %>%
      mutate(
        dyad = paste0(deployid, "_", contact_id),
        distance = distance / 1000
      ) %>%
      dplyr::select(-start_deployment_deploy, -end_deployment_deploy, 
                    -start_deployment_contact, -end_deployment_contact)
    
    # Fit zero-inflated gamma model
    zigamma_contact_mod <- glmmTMB(
      formula = contact_count ~ distance + (1 | dyad),
      ziformula = ~distance + (1 | dyad),
      data = complete_contact_summary_Age_Sex,
      family = ziGamma(link = "log")
    )
    
    return(zigamma_contact_mod)
  }
  
  if(distance_interval == "daily") {
    
    ###Merge full grid with observed contacts
    complete_data <- dyad_date_grid %>%
      left_join(daily_contact_summary, 
                by = c("deployid", "contact_id", "date")) %>%
      mutate(
        daily_contact = replace_na(daily_contact, 0),
        contact_count = replace_na(contact_count, 0)
      )
    
    
    # Filter to only dates within deployment period
    
    complete_contact_summary <- complete_data %>% 
      merge(., collar_metadata, by.x = "deployid", by.y = "animal_id") %>%
      dplyr::select(-Sex, -Age_Class) %>%
      dplyr::rename(start_deployment_deploy = start_deployment, 
                    end_deployment_deploy = end_deployment) %>%
      merge(., collar_metadata, by.x = "contact_id", by.y = "animal_id") %>%
      dplyr::select(-Sex, -Age_Class) %>%
      dplyr::rename(start_deployment_contact = start_deployment, 
                    end_deployment_contact = end_deployment) %>%
      filter(date >= start_deployment_deploy &
               (is.na(end_deployment_deploy) | date <= end_deployment_deploy),
             date >= start_deployment_contact &
               (is.na(end_deployment_contact) | 
                  date <= end_deployment_contact)) %>%
      mutate(
        dyad = paste0(deployid, "_", contact_id),
        distance = distance / 1000
      ) %>%
      dplyr::select(-start_deployment_deploy, -end_deployment_deploy, 
                    -start_deployment_contact, -end_deployment_contact)
    
    
    # Fit zero-inflated gamma model
    
    zigamma_contact_mod <- glmmTMB(
      formula = contact_count ~ distance + (1 | dyad),
      ziformula = ~distance + (1 | dyad),
      data = complete_contact_summary,
      family = ziGamma(link = "log")
    )
    
    return(zigamma_contact_mod)
  }
}

####
# create_study_area_grid: A function to create metapopulation grid across a 
# study area
create_study_area_grid <- function(median_mvmt, #locations, 
                                   buffer, 
                                   common_crs, 
                                   mask_file = NULL, 
                                   square = FALSE, 
                                   RDS = FALSE, 
                                   SHP = FALSE) {
  
  require(tidyverse)
  require(terra)
  require(sf)
  
  #Create study area extent and grid
  pix.size <- round(median_mvmt[which(median_mvmt$Sex == "Female" & 
                                        median_mvmt$Age_Class == "Adult"), 
                                "est"],0)
  #Load mask from RDS
  if(!is.null(mask_file) & RDS == TRUE)
  {boundary <- readRDS(mask_file)}
  
  #Load mask from SHP
  if(!is.null(mask_file) & SHP == TRUE)
  {boundary <- terra::vect(mask_file)}
  
  #Ensure projection of mask
  terra::crs(boundary) <- common_crs
  
  grid_cells <- sf::st_as_sf(boundary, coords = c("x","y")) %>% 
    sf::st_set_crs(common_crs) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    st_buffer(buffer) %>% 
    st_make_grid(cellsize = pix.size, square = square)
  
  grid_cells <- st_sf(geometry = grid_cells)
  grid_cells <- vect(grid_cells)
  
  #Add centroids
  centroids <- centroids(grid_cells)
  
  # Extract the X and Y coordinates of the centroids
  centroid_coords <- geom(centroids)[, c("x", "y")]
  
  # Add the centroid coordinates as new columns in the grid_cells SpatVector
  grid_cells$centroid_x <- centroid_coords[, "x"]
  grid_cells$centroid_y <- centroid_coords[, "y"]
  
  final_grid <- grid_cells
  
  #crop grid cells by boundary
  final_grid <- terra::crop(grid_cells, boundary) %>% 
    tidyterra::mutate(grid_id = 1:nrow(.))
  
  return(final_grid)
  
}

####
# pixel_options: A function to create a list of destination pixels
pixel_options <- function(grid,
                          movement_params, 
                          age_sex = "Adult male", 
                          threshold = 0.99, 
                          centroid_decay_model = NULL){
  
  require(lme4)
  require(dplyr)
  
  #Create pairwise distance df for selecting next pixel to move to
  # Assuming grid is already loaded as a SpatVector
  # Convert to an sf object if necessary
  # Remove the original polygon geometry column if necessary
  centroids_sf <- st_set_geometry(st_as_sf(grid), "geometry")
  centroids_sf_df <- st_set_geometry(centroids_sf, NULL)
  
  # Convert centroid_x and centroid_y into point geometries
  centroids_sf <- st_as_sf(
    centroids_sf_df,
    coords = c("centroid_x", "centroid_y"),
    crs = st_crs(centroids_sf)  # Use the same CRS as the original object
  )
  
  # Extract coordinates and density probabilities from the centroids
  centroids_df <- centroids_sf %>%
    st_coordinates() %>%
    as.data.frame() %>%
    mutate(pixel_id = as.numeric(rownames(centroids_sf))) %>% 
    dplyr::select(., X, Y)
  
  colnames(centroids_df) <- c("X", "Y")
  
  # Calculate pairwise distances using dist function (in meters)
  dist_matrix <- as.matrix(dist(centroids_df))
  
  # Convert the distance matrix to a data frame
  dist_df <- as.data.frame(as.table(dist_matrix))
  
  # Rename columns
  colnames(dist_df) <- c("pixel_1", "pixel_2", "distance")
  
  # Add density columns for both pixels
  dist_df <- dist_df %>%
    mutate(
      pixel_1 = as.numeric(as.character(pixel_1)),
      pixel_2 = as.numeric(as.character(pixel_2))
    )
  
  dist_df <- dist_df %>%
    mutate(dist_HR_centroid_prob_M = predict(
      centroid_decay_model$male_centroid_distance_decay,
      newdata = data.frame(distance_to_centroid = distance/1000),
      type = "response"
    ),
    dist_HR_centroid_prob_F = predict(
      centroid_decay_model$female_centroid_distance_decay,
      newdata = data.frame(distance_to_centroid = distance/1000), 
      type = "response")
    )
  return(dist_df)
}

####
# add_raster_to_grid: A function to add pixel values for a particular variable 
# derived from a raster layer.

add_raster_to_grid <- function(grid, raster_to_add, field_to_add, 
                               common_crs = "EPSG:32615", display_plot = TRUE){
  require(sf)
  require(terra)
  require(dplyr)
  
  #Bring in raster
  raster <- terra::rast(raster_to_add)
  
  #Match projection with grid_cells object
  raster <- terra::project(x = raster, common_crs)
  
  # Extract and summarize raster values (mean) within each polygon
  extracted_values <- terra::extract(raster, grid, fun = mean, na.rm = TRUE)
  colnames(extracted_values)[2] <- field_to_add
  
  extracted_values[is.nan(extracted_values[,2]),2] <- 0
  # Add the extracted values as a new column to the sf object
  grid[,field_to_add] <- extracted_values[,field_to_add]
  
  if(display_plot == TRUE) {terra::plot(grid, y = field_to_add, 
                                        type = "continuous")}
  return(grid)
}

####
# add_centroid_count_to_grid: A function to calculate number of point locations 
# within grid cells
add_centroid_count_to_grid <- function(shapefile_path, 
                                       name, 
                                       grid, 
                                       crs = "EPSG:32615", 
                                       display_plot = TRUE, 
                                       prefiltered = TRUE){
  
  # Step 1: Read in the shapefile
  shp <- vect(shapefile_path, crs = crs)  # Using terra to load shapefile
  
  # Step 2: Filter polygons by an attribute (if needed)
  if(prefiltered == FALSE)
  {attribute_field_1 <- "HMSTD_CD1"  # Replace with your field name
  text_string_1 <- "H"  # Replace with your target text
  attribute_field_2 <- "PR_TYP_NM1"  # Replace with your field name
  text_string_2 <- "RESIDENTIAL"  # Replace with your target text
  filtered_shp <- shp[shp[[attribute_field_1]] == text_string_1 & 
                        shp[[attribute_field_2]] == text_string_2, ]}
  
  if(prefiltered == TRUE) {filtered_shp <- shp}
  
  # Step 3: Generate centroids for the filtered polygons
  centroids <- centroids(filtered_shp)
  
  # Step 4: Count centroids in each grid cell
  # Check the intersection of centroids with grid cells
  intersections <- terra::extract(grid, centroids)
  
  # Count the number of centroids in each polygon
  centroid_counts <- intersections %>%
    filter(!is.na(grid_id)) %>% 
    group_by(grid_id) %>%
    dplyr::summarize(count = n(), .groups = "drop")
  
  # Add the count to the EC.grid SpatVector
  grid$centroid_count <- 0  # Initialize with 0
  grid$centroid_count[centroid_counts$grid_id] <- centroid_counts$count
  
  names(grid)[names(grid) == "centroid_count"] <- name
  
  if(display_plot == TRUE) {terra::plot(grid, y = name, 
                                        col = terrain.colors(10), 
                                        type = "continuous")}
  return(grid)
}

####
# set_simulation_parameters: A collector function to define all transmission 
# parameters for a disease simulation

set_simulation_parameters <- function(n.sims,
                                      init.prev, 
                                      init.I, 
                                      init.group, 
                                      vary_epi_params = FALSE)
{
  require(whitetailedSIRS)
  
  if(vary_epi_params == TRUE){
    #Read in elicitation data (Rosenblatt et al. 2024)
    elicitation_data <- whitetailedSIRS::draw_elicitation_samples(
      nsamples = n.sims)
    
    #Create tibble with all parameters used
    sim_params <- tibble(sim = 1:n.sims,
                         init.prev = init.prev,
                         init.I = init.I,
                         init.group = init.group,
                         
                         alpha_immunity = 1 / get_EE_param_vals(data = 
                            elicitation_data, my_param = 'Temporary Immunity'),
                         
                         C_nu_deer = 10^5.6 * get_EE_param_vals(data = 
                            elicitation_data, my_param = "Viral Load"),
                         
                         r_deer = get_EE_param_vals(data = elicitation_data,
                            my_param = "Dose-Response"), 
                         
                         t_contact_deer_deer = get_EE_param_vals(data = 
                            elicitation_data, 
                            my_param = "Deer Proximity Duration (minutes)"),
                         
                         sigma_aero_deer_deer = calc_sigma_aero(
                           C_nu = C_nu_deer,
                           t_contact = t_contact_deer_deer / 60, 
                           r = r_deer),
                         
                         epsilon_dc = get_EE_param_vals(data = elicitation_data, 
                           my_param = "Direct Contact Probability"),
                         
                         sigma_dc_deer_deer = calc_sigma_dc(C_nu = C_nu_deer, 
                                                            nsamples = n.sims),
                         
                         C_nu_human = rnorm(n = n.sims, 
                                            mean = 10^5.6, 
                                            sd = 10^1.2),
                         
                         t_contact_deer_human = get_EE_param_vals(data = 
                           elicitation_data, 
                           my_param = 
                           "Deer-Human Proximity Duration, Suburban (minutes)"),
                         
                         #Estimate duration of 
                         #human-deer proximity event in suburban context...
                         sigma_aero_deer_human = calc_sigma_aero(ER = 0.53, 
                          C_nu = C_nu_human, 
                          t_contact = t_contact_deer_human / 60, 
                          r = r_deer, 
                          nsamples = n.sims)#... and calculate infection 
                         #probability given the duration of a human-deer proximity
                         #event.
    )
    
    if(sum(sim_params$init.prev > 0 & sim_params$init.I >0) >=1 |
       sum(sim_params$init.prev > 0 & sim_params$init.group >0) >=1 |
       sum(sim_params$init.I > 0 & sim_params$init.group >0) >=1 ) 
    {stop(print("There is one or more parameter sets that have non-zero values
                  for both init.prev, init.I, and init.group. Adjust these 
                arguments to ensure that each simulation uses init.prev OR 
                init.I OR init.group, but not combinations of the three"))}
  }
  
  if(vary_epi_params == FALSE){
    #Create tibble with all parameters used
    sim_params <- tibble(sim = 1:n.sims,
                         init.prev = init.prev,
                         init.I = init.I,
                         init.group = init.group)
  }
  
  return(sim_params)
}

#### 
# run_simulation_local_density_est_fast: A function to run one simulation of 
# disease transmission. Output focuses on a summary, rather than large tibbles 
# of infected individuals
run_simulation_local_density_est_fast <- function(a = i,
                                                  starting.values,
                                                  sim.input,
                                                  grid,
                                                  centroids_df,
                                                  mvmt_dist_shape,
                                                  cand_pixels,
                                                  contact_distance_df,
                                                  contact_distance_sigma,
                                                  sigma_infection_prob,
                                                  median_mvmt,
                                                  stop_I_0,
                                                  prop_residence_contact,
                                                  rate_residence_contact) {
  
  # --- Initialize population ---
  individuals_t <- initiate_population(
    n_pop = sim.input$pop_size, 
    avg_group_size = sim.input$ave_group_size, 
    sex_ratio = sim.input$sex_ratio,
    grid = grid, roster = sim.input$collar_roster, 
    origin_cell_pref = sim.input$origin_probs, 
    prop.infected = starting.values$init.prev[a], 
    number_init_infected = starting.values$init.I[a], 
    number_group_infected = starting.values$init.group[a], 
    centroids_df = centroids_df
  )
  
  n_steps <- 91
  
  infection_summary_list <- vector("list", n_steps)
  
  start_cell <- unique(individuals_t$Cell[individuals_t$I == 1])

  cell_use_list    <- vector("list", n_steps)
  infected_group_list  <- vector("list", n_steps)
  compartment_list <- vector("list", n_steps)
  group_size_list <- vector("list", n_steps)
  
  # --- Initial state ---
  cell_use_list[[1]] <- individuals_t[, c("Group","Cell")]
  cell_use_list[[1]]$time <- 0
  
  group_zero <- individuals_t[individuals_t$I == 1, ] |>
    dplyr::select(Group) |>
    dplyr::mutate(time = 0)
  
  group_size_list[[1]] <- individuals_t[individuals_t$S == 1, ] |>
    dplyr::count(Group, name = "N") |>
    dplyr::mutate(time = 0)
  
  compartment_list[[1]] <- data.frame(
    S = mean(individuals_t$S),
    I = mean(individuals_t$I),
    R = mean(individuals_t$R),
    time = 0,
    sim = a
  )
  
  # --- FAST infection duration tracking (ID SAFE) ---
  infection_days <- setNames(
    integer(nrow(individuals_t)),
    individuals_t$ID
  )
  
  infection_days[individuals_t$ID[individuals_t$I == 1]] <- 0
  
  # helper function (lightweight infection summary)
  
  summarise_infections <- function(individuals_t, 
                                   idx, 
                                   time, 
                                   sim_id, 
                                   start_cell) {
    
    if (!any(idx)) {
      return(list(
        cases = 0,
        groups = 0,
        time = time,
        sim = sim_id,
        origin_cell = start_cell
      ))
    }
    
    list(
      cases = sum(idx),
      groups = length(unique(individuals_t$Group[idx])),
      time = time,
      sim = sim_id,
      origin_cell = start_cell
    )
  }
  
  # time step loop
  for (i in 1:n_steps) {
    
    # 1. Move groups
    individuals_t <- move_groups(
      individuals_df = individuals_t, 
      grid = grid,
      age_sex_params = mvmt_dist_shape, 
      candidate_cells = cand_pixels, 
      median_mvmt = median_mvmt,
      pixel_selection_df = sim.input$selection_probs,
      distance_from_centroid = sim.input$centroid_decay,
      centroids_df = centroids_df
    )
    
    # 2. Cell tracking
    cell_use_list[[i]] <- cbind(
      individuals_t[, c("Group","Cell")],
      time = i
    )
    
    # 3. Group size tracking
    group_size_list[[i]] <- individuals_t[individuals_t$S == 1, ] |>
      dplyr::count(Group, name = "N") |>
      dplyr::mutate(time = i)
    
    # 4. FOI
    individuals_t <- FOI_builder(
      candidate_pixels = cand_pixels,
      sigma_infection_prob = sigma_infection_prob,
      grid = grid,
      contact_distance_df = contact_distance_df,
      contact_distance_sigma = contact_distance_sigma,
      individuals_df = individuals_t,
      starting.values = starting.values[a,],
      prop_residence_contact = prop_residence_contact,
      rate_residence_contact = rate_residence_contact
    )
    
    # infection summary (NO DATAFRAMES CREATED)
    idx <- ((individuals_t$new_infect_from_deer == 1 |
               individuals_t$new_infect_from_humans == 1) &
              individuals_t$I == 1)
    
    infection_summary_list[[i]] <- summarise_infections(
      individuals_t,
      idx,
      i,
      a,
      start_cell
    )
    
    infected_group_list[[i]] <- individuals_t[individuals_t$I == 1,] |>
      dplyr::select(Group) |>
      dplyr::mutate(time = i)
    
    current_ids <- individuals_t$ID
    
    missing_ids <- setdiff(current_ids, names(infection_days))
    if (length(missing_ids) > 0) infection_days[missing_ids] <- 0
    
    infection_days[current_ids[individuals_t$I == 1]] <-
      infection_days[current_ids[individuals_t$I == 1]] + 1
    
    infection_days[current_ids[individuals_t$I == 0]] <- 0
    
    recover_ids <- names(infection_days)[infection_days == 7]
    
    individuals_to_recover <- data.frame(ID = recover_ids)
    
    # 5. State transitions
    individuals_t <- perform_state_transitions(
      individuals_t, 
      starting.values$alpha_immunity[a], 
      force_recovery = TRUE, 
      no_immune_loss = TRUE,
      recovering_individuals = individuals_to_recover
    )
    
    # infection summary (NO DATAFRAMES CREATED)
    idx <- ((individuals_t$new_infect_from_deer == 1 |
               individuals_t$new_infect_from_humans == 1) &
              individuals_t$I == 1)
    
    infection_summary_list[[i]] <- summarise_infections(
      individuals_t,
      idx,
      i,
      a,
      start_cell
    )
    
    # 6. Trim object
    individuals_t <- individuals_t[, c("Cell", "ID", "Sex", "clone_of", "Group",
                                       "Centroid_X", "Centroid_Y", 
                                       "Centroid_cell",
                                       "Current_X", "Current_Y", "prev_cell",
                                       "S", "I", "R", "FOI", "FOI_human")]
    
    # 7. Compartment sizes
    compartment_list[[i]] <- data.frame(
      S = mean(individuals_t$S),
      I = mean(individuals_t$I),
      R = mean(individuals_t$R),
      time = i,
      sim = a
    )
    
    # Early stop
    if (stop_I_0) {
      if (compartment_list[[i]]$I == 0 | compartment_list[[i]]$S == 0) break
      }
    }
  
  # FINAL COMBINES

  cell_use_tracker  <- dplyr::bind_rows(cell_use_list)
  group_sizes       <- dplyr::bind_rows(group_size_list)
  infected_groups   <- dplyr::bind_rows(infected_group_list) %>% 
    dplyr::select(Group) %>% unique()
  compartment_sizes <- dplyr::bind_rows(compartment_list)
  infection_summary <- dplyr::bind_rows(infection_summary_list)
  
  # WEEKLY SUMMARY (UPDATED INPUT)
  weekly_group_size <- group_sizes %>%
    mutate(week = floor(time / 8)+1)
  
  weekly_infections <- infection_summary %>%
    mutate(week = floor(time / 8)+1)
  
  infected_cells <- cell_use_tracker %>%
    inner_join(infected_groups, by = c("Group")) %>%
    distinct(Cell)
  
  weekly_results <- cell_use_tracker %>%
    inner_join(infected_cells, by = "Cell") %>%
    distinct(Group) %>% 
    left_join(weekly_group_size) %>%
    group_by(time) %>%
    summarise(local_N = sum(N, na.rm = TRUE), .groups = "drop") %>%
    mutate(week = floor(time / 8)+1) %>% 
    group_by(week) %>% 
    reframe(local_N = first(local_N)) %>% 
    left_join(., 
      weekly_infections %>%
        group_by(week) %>%
        summarise(new_cases = sum(cases), .groups = "drop"),
      by = "week"
    ) %>%
    mutate(
      new_cases = ifelse(is.na(new_cases), 0, new_cases),
      rate = new_cases / local_N
    )

  if(nrow(infected_groups) == 1) {spatial_spread = tibble(week = 1, spatial_spread = 0)}
  if(nrow(infected_groups) > 1)  {spatial_spread <- dplyr::bind_rows(
    infected_group_list) %>%
    unique() %>% #removes duplicates within group
    mutate(week = floor(time / 8)+1) %>%
    filter(!Group %in% group_zero$Group) %>% 
    inner_join(unique(cell_use_tracker), by = c("Group", "time")) %>% 
    mutate(origin_pixel = unique(infection_summary$origin_cell)) %>% 
    unique() %>% 
    inner_join(cand_pixels %>% dplyr::select(pixel_1, pixel_2, distance), 
               by = c("origin_pixel" = "pixel_1", "Cell" = "pixel_2")) %>% 
    group_by(week) %>% 
    reframe(spatial_spread = max(distance/1000)) %>% 
    arrange(week) %>% 
    mutate(spatial_spread = pmax(spatial_spread - dplyr::lag(spatial_spread, 
                                                             default = 0),0))}
    

  
  output <- tibble(simulation_id = a, 
                   total_cases = sum(infection_summary$cases),
                   end = last(compartment_sizes$time),
                   compartment_sizes %>% reframe(time_max_prev = 
                                                   time[which.max(I)]),
                   origin_pixel = unique(infection_summary$origin_cell),
                   n_groups = nrow(infected_groups),
                   max_daily_incidence = max(infection_summary$cases),
                   incidence_category = case_when(
                     sum(infection_summary$cases) == 0 ~ "Failed",
                     sum(infection_summary$cases) > 0 & 
                       max(infection_summary$cases) == 1 ~ "Stutter",
                     max(infection_summary$cases) > 0 ~ "Outbreak"
                   ),
                   max_local_epi_rate = max(weekly_results$rate, na.rm = T),
                   max_spatial_spread = max(spatial_spread$spatial_spread),
                   end_seroprev = last(compartment_sizes$R)
  )
  
  # RETURN (LIGHTWEIGHT)
  list(
    infection_summary = infection_summary,
    compartment_sizes = compartment_sizes,
    local_epi_growth_rate_weekly = weekly_results,
    output = output
  )
}

####
# run_all_simulations_local_density_est: A wrapper function to run multiple 
# iterations of a disease simulation

run_all_simulations_local_density_est <- function(starting.values,
                                sim.input,
                                grid,
                                centroids_df,
                                mvmt_dist_shape,
                                cand_pixels,
                                contact_distance_df,
                                contact_distance_sigma,
                                sigma_infection_prob,
                                median_mvmt,
                                stop_I_0 = TRUE,
                                prop_residence_contact,
                                rate_residence_contact,
                                num_cores = 2) { # Added num_cores parameter
  
  # Use mclapply for parallel processing of simulations
  # The 'mc.cores' argument specifies the number of cores to use.
  
  results <- mclapply(1:nrow(starting.values), function(i) {
    message("Running simulation ", i, " of ", nrow(starting.values))
    
    # Run the simulation and handle errors
    result <- tryCatch(
      run_simulation_local_density_est_fast(
        a = i,
        starting.values = starting.values,
        sigma_infection_prob = sigma_infection_prob,
        sim.input = sim.input,
        grid = grid,
        centroids_df = centroids_df,
        mvmt_dist_shape = mvmt_dist_shape,
        cand_pixels = cand_pixels,
        contact_distance_df = contact_distance_df,
        contact_distance_sigma = contact_distance_sigma,
        median_mvmt = median_mvmt,
        stop_I_0 = stop_I_0,
        prop_residence_contact = prop_residence_contact,
        rate_residence_contact = rate_residence_contact
      ),
      error = function(e) {
        message("Error in simulation ", i, ": ", e$message)
        list(error = e$message)  # Return error details
      }
    )
    
    return(result)
  }, mc.cores = num_cores, mc.preschedule = FALSE) # Pass num_cores to mclapply
  
  # Combine successful results into a single list of data frames
  combined_results <- list(
    infection_summary = bind_rows(lapply(results,
                            function(res) if (!is.null(res$infection_summary)) 
                              res$infection_summary),
                               .id = "simulation_id"),
    compartment_sizes = bind_rows(lapply(results,
                            function(res) if (!is.null(res$compartment_sizes)) 
                              res$compartment_sizes),
                                  .id = "simulation_id"),
    local_epi_growth_rate_weekly = bind_rows(lapply(results,
                  function(res) if (!is.null(res$local_epi_growth_rate_weekly)) 
                    res$local_epi_growth_rate_weekly),
                                              .id = "simulation_id"),
    output = bind_rows(lapply(results, function(res) if (!is.null(res$output)) 
      res$output),
                                             .id = "simulation_id"))
  
  # Return combined results
  return(combined_results)
}

####
# initiate_population: Generate initial deer population for simulation
initiate_population <- function(
    n_pop,
    avg_group_size,
    sex_ratio,
    grid,
    roster,
    origin_cell_pref,
    prop.infected = 0,
    number_init_infected = 0,
    number_group_infected = 0,
    centroids_df) 
  
{
  
  # Consistency check for infection settings
  if (sum(c(prop.infected > 0, 
            number_init_infected > 0, 
            number_group_infected > 0) > 1)) {
    stop("Specify only one infection parameter: 
         proportion OR number of individuals OR number of groups.")
  }
  
  # Calculate sex split
  num_males <- round(n_pop * sex_ratio / (1 + sex_ratio))
  num_females <- n_pop - num_males
  
  # Fast group size generator: vectorized
  #Generate group sizes for males and females using truncated 
  # Poisson distribution
  generate_group_sizes <- function(total_count, avg_size) {
    group_sizes <- c()
    while (sum(group_sizes) < total_count) {
      group_sizes <- c(group_sizes, rtpois(1, lambda = avg_size, a = 0))
    }
    if (sum(group_sizes) > total_count) {
      group_sizes[length(group_sizes)] <- total_count - 
        sum(group_sizes[-length(group_sizes)])
    }
    return(group_sizes)
  }
  
  male_group_sizes <- generate_group_sizes(num_males, avg_group_size)
  female_group_sizes <- generate_group_sizes(num_females, avg_group_size)
  
  # Assign IDs & groups, vectorized
  male_individuals <- data.frame(
    ID = seq_len(num_males),
    Sex = "Male",
    Group = rep(seq_along(male_group_sizes), times = male_group_sizes)
  )
  
  female_individuals <- data.frame(
    ID = (num_males + 1):n_pop,
    Sex = "Female",
    Group = rep(seq_along(female_group_sizes) + length(male_group_sizes),
                times = female_group_sizes)
  )
  
  individuals <- dplyr::bind_rows(male_individuals, female_individuals)
  
  # Assign what collared animal each simulated group is based upon
  unique_female_groups <- individuals %>%
    filter(Sex == "Female") %>% 
    distinct(Group)
  
  unique_male_groups <- individuals %>%
    filter(Sex == "Male") %>% 
    distinct(Group)

  unique_female_groups <- unique_female_groups %>%
    mutate(
      clone_of = sample(
        roster$animal_id[which(roster$Sex == "F")],
        size = n(),
        replace = TRUE
      )
    )
  
  unique_male_groups <- unique_male_groups %>%
    mutate(
      clone_of = sample(
        roster$animal_id[which(roster$Sex == "M")],
        size = n(),
        replace = TRUE
      )
    )
  
  unique_groups <- rbind(unique_female_groups, unique_male_groups)
  
  # Assign centroid cell to each group based off of the collared animal 
  # the group's characteristics are based on.
  group_grid <- unique_groups %>%
    rowwise() %>%
    mutate(
      grid_id = sample(
        origin_cell_pref$grid_id,
        size = 1,
        replace = TRUE,
        prob = origin_cell_pref[[clone_of]]
      )
    ) %>%
    ungroup() %>%
    left_join(centroids_df[c("grid_id", "X", "Y")], by = "grid_id")
  
  
  # Merge back, add coordinates & states
  individuals <- individuals %>%
    left_join(group_grid, by = "Group") %>%
    rename(
      Centroid_cell = grid_id,
      Centroid_X = X,
      Centroid_Y = Y
    ) %>%
    mutate(
      Cell = Centroid_cell,
      Current_X = Centroid_X,
      Current_Y = Centroid_Y,
      prev_cell = NA,
      S = 1, I = 0, R = 0
    )
  
  # Infection assignment — vectorized
  if (prop.infected > 0) {
    infected_indices <- sample(
      seq_len(nrow(individuals)),
      size = round(prop.infected * nrow(individuals)),
      replace = FALSE
    )
    individuals$I[infected_indices] <- 1
    individuals$S[infected_indices] <- 0
  }
  
  if (number_init_infected > 0) {
    infected_indices <- sample(seq_len(nrow(individuals)), 
                               size = number_init_infected)
    individuals$I[infected_indices] <- 1
    individuals$S[infected_indices] <- 0
  }
  
  if (number_group_infected > 0) {
    infected_groups <- sample(unique(individuals$Group), 
                              size = number_group_infected)
    individuals <- individuals %>%
      mutate(
        I = if_else(Group %in% infected_groups, 1L, I),
        S = if_else(I == 1, 0L, S)
      )
  }
  
  return(individuals)
}



####
# move_groups: A function that initiates group movement for each time step
move_groups <- function(
    individuals_df, 
    grid, 
    age_sex_params, 
    candidate_cells, 
    median_mvmt, 
    pixel_selection_df, 
    distance_from_centroid,
    centroids_df
) {
  
  # 1. Precompute exit pixel distance ONCE
  exit_pixel_dist <- median_mvmt %>%
    filter(Sex %in% c("F","Female"), Age_Class == "Adult") %>%
    pull(est) %>%
    round(0)
  
  # 2. Compute unique group centroids
  group_centroids_df <- individuals_df %>%
    group_by(Group) %>%
    summarise(
      Sex = first(Sex),
      clone_of = first(clone_of),
      Centroid_X = first(Centroid_X),
      Centroid_Y = first(Centroid_Y),
      Centroid_cell = first(Centroid_cell),
      Current_X = first(Current_X),
      Current_Y = first(Current_Y),
      Cell = first(Cell),
      prev_cell = first(prev_cell),
      .groups = "drop"
    )
  
  # 3. Vectorized random distance sampling
  updated_group_df <- group_centroids_df %>%
    mutate(
      age_sex_key = if_else(Sex == "Male", "Adult male", "Adult female")
    ) %>%
    left_join(age_sex_params, by = c("age_sex_key" = "Age_Sex")) %>%
    mutate(
      random_distance = rgamma(n(), shape = shape, rate = rate)
    )
  
  # Check for NA in random_distance quickly
  if (anyNA(updated_group_df$random_distance)) {
    stop("Random distance generation resulted in NA.")
  }
  
  # 4. Split groups: stay vs move
  stay <- updated_group_df %>%
    filter(random_distance <= exit_pixel_dist) %>%
    mutate(prev_cell = Cell) %>% 
    dplyr::select(Group, Sex, clone_of, Centroid_X, Centroid_Y, Centroid_cell, 
                  Current_X, Current_Y, Cell, prev_cell)
  
  move <- updated_group_df %>%
    filter(random_distance > exit_pixel_dist)
  
  # Edge case check
  if (nrow(stay) == 0 | nrow(stay) == nrow(updated_group_df)) {
    stop("All groups stayed put or moved. Check random distance generation.")
  }
  
  # 5. Candidate options for moving groups. is_normal indicates that the random 
  # distance is within a grid cell diameter from the true distance
  options <- move %>%
    left_join(candidate_cells, by = c("Cell" = "pixel_1"), 
              relationship = "many-to-many") %>%
    mutate(
      is_normal = abs(distance - random_distance) <= exit_pixel_dist
    )
  
  options_normal <- options %>% filter(is_normal)
  
  # Options_exceed indicate animals making large movements, that would take 
  # them out of the study area. These groups are forced to take the next best 
  # option
  options_exceed <- options %>%
    filter(!(Group %in% options_normal$Group)) %>%
    group_by(Group) %>%
    filter(if (n() == 0) FALSE else 
      distance >= (max(distance, na.rm = TRUE) - exit_pixel_dist)) %>%
    ungroup()
  
  combined_options <- bind_rows(options_normal, options_exceed)
  
  rsf_lookup <- pixel_selection_df %>%
    pivot_longer(-grid_id, names_to = "clone_of", values_to = "RSF_prob")
  
  candidates <- combined_options %>%
    left_join(rsf_lookup, by = c("pixel_2" = "grid_id", "clone_of")) %>%
    mutate(
      dist_HR_centroid_prob = if_else(Sex == "Male", dist_HR_centroid_prob_M, 
                                      dist_HR_centroid_prob_F)
    ) %>%
    dplyr::select(-dist_HR_centroid_prob_M, -dist_HR_centroid_prob_F)
  
  # 6. Conspecific density effect
  cell_counts <- individuals_df %>%
    dplyr::count(Cell, name = "n_individuals_in_cell")
  
  density_effect <- cell_counts %>%
    mutate(
      pixel_2 = Cell,
      density_prob = 1 / (n_individuals_in_cell + 1)
    ) %>%
    dplyr::select(pixel_2, density_prob)
  
  candidates <- candidates %>%
    left_join(density_effect, by = "pixel_2") %>%
    mutate(density_prob = coalesce(density_prob, 1))  # 1 = lowest density
  
  # 7. Final combined probability
  candidates <- candidates %>%
    mutate(product_prob = RSF_prob * dist_HR_centroid_prob * density_prob)
  
  # 8. Sample destination per group
  moved <- candidates %>%
    group_by(Group) %>%
    slice_sample(n = 1, weight_by = product_prob) %>%
    ungroup() %>%
    mutate(prev_cell = Cell) %>% 
    mutate(Cell = pixel_2) %>%
    dplyr::select(Group, Sex, clone_of, Centroid_X, Centroid_Y, Centroid_cell, 
                  Current_X, Current_Y, Cell, prev_cell)
  
  # 9. Add new XY by joining precomputed centroids
  moved <- moved %>%
    left_join(centroids_df[,c("X", "Y", "grid_id")], 
              by = c("Cell" = "grid_id")) %>%
    mutate(Current_X = X, Current_Y = Y) %>%
    dplyr::select(-X, -Y)
  
  # 10. Combine stay + moved
  updated_group_centroids <- bind_rows(stay, moved)
  
  if (anyNA(updated_group_centroids$Cell)) {
    stop("NA in new pixels selected.")
  }
  
  # 11. Update full individuals_df
  individuals_df_2 <- individuals_df %>%
    dplyr::select(-Centroid_X, -Centroid_Y, -Centroid_cell, -Current_X, 
                  -Current_Y, -Cell, -prev_cell) %>%
    left_join(updated_group_centroids, by = c("Group", "Sex", "clone_of")) %>%
    dplyr::select(ID, Sex, clone_of, Group, Centroid_X, Centroid_Y, 
                  Centroid_cell, Current_X, Current_Y, Cell, prev_cell, S, I, R)
  
  if (anyNA(individuals_df_2$Cell)) {
    stop("An individual moved to an NA cell.")
  }
  
  return(individuals_df_2)
}

####
# FOI_builder: A function to draw values for estimating FOI
# using a single value (0.018)
FOI_builder <- function(candidate_pixels, 
                                     sigma_infection_prob,
                                     grid, 
                                     contact_distance_df,
                                     contact_distance_sigma,
                                     individuals_df, 
                                     starting.values, 
                                     prop_residence_contact = 0, 
                                     rate_residence_contact = 0) {
  require(dplyr)
  require(whitetailedSIRS)
  
  # Aggregate individuals by cell, ignoring role of sex as contact rates and 
  # susceptibility are the same
  infected_locations <- individuals_df %>% 
    filter(I > 0) %>% 
    dplyr::select(ID, Cell)
  
  susceptible_df <- candidate_pixels[, c("pixel_1","pixel_2", "distance")] %>%
    filter(pixel_1 %in% unique(individuals_df$Cell[which(individuals_df$S==1)]), 
           pixel_2 %in% unique(infected_locations$Cell)) %>% 
    left_join(., individuals_df[which(individuals_df$S == 1),], 
              by = c("pixel_1" = "Cell"), 
              relationship = "many-to-many") %>% 
    dplyr::select(ID, pixel_1, pixel_2, distance) %>% 
    rename("ID_S" = "ID") %>% 
    left_join(., infected_locations, 
              by = c("pixel_2" = "Cell"), 
              relationship = "many-to-many") %>% 
    rename("ID_I" = "ID") %>% 
    mutate(distance_km = round(distance/1000, digits = 3)) %>% 
    dplyr::select(ID_S, ID_I, distance_km) %>% 
    left_join(., contact_distance_df, by = "distance_km") %>% 
    mutate(contact = case_when(is.na(contact) ~ 0, TRUE ~ contact)) %>% 
    mutate(proximity_rate_pred = contact,
           proximity_shape = 1 / (contact_distance_sigma^2),
           proximity_scale = proximity_rate_pred / proximity_shape,
           proximity_rate_draw = rgamma(n(), shape = proximity_shape, 
                                        scale = proximity_scale),
           beta_pairwise = proximity_rate_draw * sigma_infection_prob
    ) %>% 
    dplyr::select(ID_S, ID_I, distance_km, beta_pairwise) %>% 
    group_by(ID_S) %>% 
    reframe(FOI = sum(beta_pairwise))
  
  
  individuals_df <- individuals_df %>%
    left_join(susceptible_df, by = c("ID" = "ID_S")) %>%
    mutate(FOI = case_when(is.na(FOI) ~ 0, TRUE ~ FOI))
  
  #Integrate human FOI (to be used when wanted) This needs to be fixed to be 
  #usable on an interaction by interaction basis*******
  
  individuals_df_w_FOI <- individuals_df %>% 
    merge(., grid, by.x = "Cell", by.y = "grid_id") %>%
    dplyr::select(., -"centroid_x", -"centroid_y") %>% 
    mutate(proximity_human = prop_residence_contact * 
             rate_residence_contact * Residence_count) %>% 
    mutate(sigma_aero_deer_human = 0,
           Beta_humans_aero = proximity_human * sigma_aero_deer_human,
           FOI_human = Beta_humans_aero*0.05) %>% 
    dplyr::select(-proximity_human, -sigma_aero_deer_human, -Beta_humans_aero, 
                  -Residence_count)
  
  return(individuals_df_w_FOI)
}

####
# perform_state_transitions: A function to generate state transitions for all 
# individuals in a population

perform_state_transitions <- function(population, 
                                      immunity_draws, 
                                      force_recovery = FALSE,
                                      no_immune_loss = TRUE,
                                      recovering_individuals){
  
  if(force_recovery == TRUE & no_immune_loss == FALSE){
    population <- population %>%
      mutate( # Transition events
        recover = case_when(I == 1 & 
                              ID %in% recovering_individuals$ID == 1 ~ 1, 
                            TRUE ~ 0),
        immune_loss = case_when(R == 1 & 
                                  rbinom(n(), 1, immunity_draws) == 1 ~ 1, 
                                TRUE ~ 0),
        new_infect_from_deer = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI)) == 1 ~ 1, 
                                  TRUE ~ 0),
        new_infect_from_humans = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI_human)) == 1 ~ 1, 
                                  TRUE ~ 0),
        # Update states based on transitions
        S = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 0, 
                      immune_loss == 1 ~ 1, TRUE ~ S),
        I = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 1, 
                      recover == 1 ~ 0, TRUE ~ I),
        R = case_when(recover == 1 ~ 1, immune_loss == 1 ~ 0, TRUE ~ R)
      )
  }
  
  if(force_recovery == FALSE & no_immune_loss == FALSE){
    population <- population %>%
      mutate( # Transition events
        recover = case_when(I == 1 & rbinom(n(), 1, 1 / 6) == 1 ~ 1, TRUE ~ 0),
        immune_loss = case_when(R == 1 & 
                                  rbinom(n(), 1, immunity_draws) == 1 ~ 1, 
                                TRUE ~ 0),
        new_infect_from_deer = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI)) == 1 ~ 1, 
                                  TRUE ~ 0),
        new_infect_from_humans = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI_human)) == 1 ~ 1, 
                                  TRUE ~ 0),
        
        # Update states based on transitions
        S = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 0, 
                      immune_loss == 1 ~ 1, TRUE ~ S),
        I = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 1, 
                      recover == 1 ~ 0, TRUE ~ I),
        R = case_when(recover == 1 ~ 1, immune_loss == 1 ~ 0, TRUE ~ R)
      )
  }
  
  if(force_recovery == TRUE & no_immune_loss == TRUE){
    population <- population %>%
      mutate( # Transition events
        recover = case_when(I == 1 & ID %in% 
                              recovering_individuals$ID == 1 ~ 1, TRUE ~ 0),
        new_infect_from_deer = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI)) == 1 ~ 1, 
                                  TRUE ~ 0),
        new_infect_from_humans = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI_human)) == 1 ~ 1, 
                                  TRUE ~ 0),
        
        # Update states based on transitions
        S = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 0, 
                      TRUE ~ S),
        I = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 1, 
                      recover == 1 ~ 0, TRUE ~ I),
        R = case_when(recover == 1 ~ 1, TRUE ~ R)
      )
  }
  
  if(force_recovery == FALSE & no_immune_loss == TRUE){
    population <- population %>%
      mutate( # Transition events
        recover = case_when(I == 1 & rbinom(n(), 1, 1 / 6) == 1 ~ 1, TRUE ~ 0),
        new_infect_from_deer = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI)) == 1 ~ 1, 
                                  TRUE ~ 0),
        new_infect_from_humans = case_when(S == 1 & I == 0 & R == 0 & 
                                  rbinom(n(), 1, 1 - exp(-FOI_human)) == 1 ~ 1, 
                                  TRUE ~ 0),
        
        # Update states based on transitions
        S = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 0, 
                      TRUE ~ S),
        I = case_when(new_infect_from_deer | new_infect_from_humans == 1 ~ 1, 
                      recover == 1 ~ 0, TRUE ~ I),
        R = case_when(recover == 1 ~ 1, TRUE ~ R)
      )
  }
  
  return(population)
}

####
# read_raw_GPS_fixes: A function to convert CSVs from multiple animals into a 
# list object for use with other functions.

read_raw_GPS_fixes <- function(directory, 
                               file, 
                               CSV = FALSE, 
                               SHP = TRUE, 
                               individual_files = FALSE){
  require(tidyverse)
  require(lubridate)
  require(data.table)
  
  if(CSV == TRUE & SHP == FALSE & individual_files == TRUE){
    
    #Get the list of all CSV files in the directory
    file_list <- list.files(path = directory, pattern = "*.csv", 
                            full.names = TRUE)
    
    # Function to read each CSV file
    read_csv_files <- function(file) {
      read.csv(file, stringsAsFactors = FALSE)
    }
    
    # Read all CSV files into a list of data frames
    data_list <- lapply(file_list, read_csv_files)
    
    # Renaming the list based on the value in the first row, second column of 
    # each data frame
    names(data_list) <- sapply(data_list, function(df) df[1, 2])
    
    # Define the columns to remove (e.g., "column_name1", "column_name2")
    columns_to_remove <- c("X", "age","dop","collar","collar_deployment_date", 
                           "collar_end_date","duplicate_","fast_step_", 
                           "fast_roundtrip_","nrow")
    
    # Function to remove columns from each data frame in the list
    remove_columns <- function(df, cols) {
      df[, !(names(df) %in% cols)]  # Select all columns except those in 'cols'
    }
    
    # Apply the function to each data frame in the list
    data_list <- lapply(data_list, remove_columns, cols = columns_to_remove)
    
    # Define the new column names
    new_column_names <- c("deployid","x", "y", "date_time","long","lat")
    
    # Function to rename columns
    rename_columns <- function(df, new_names) {
      colnames(df) <- new_names
      return(df)
    }
    
    # Apply the renaming function to all data frames in the list
    data_list <- lapply(data_list, rename_columns, new_names = new_column_names)
    
    # Function to calculate seconds from a reference time 
    # (e.g., the first timestamp) in each data frame
    calculate_seconds <- function(df) {
      # Convert 'timestamp' to POSIXct if it's not already
      df$date_time <- ymd_hms(df$date_time)
      
      # Calculate the seconds difference from the first timestamp 
      # (or any reference time)
      df$seconds <- as.numeric(seconds(df$date_time))
      
      return(df)
    }
    
    # Apply the function to each data frame in the list
    data_list <- lapply(data_list, calculate_seconds)
    
    # Function to filter out rows with NA values from a data frame
    remove_na_rows <- function(df) {
      df <- na.omit(df)  # Removes rows with any NA values
      return(df)
    }
    
    # Apply the function to each data frame in the list
    data_list <- lapply(data_list, remove_na_rows)
    
    # Define the new column names
    new_column_names <- c("deployid","x", "y", "date_time","longitude",
                          "latitude","t")
    
    # Function to rename columns
    rename_columns <- function(df, new_names) {
      colnames(df) <- new_names
      return(df)
    }
    
    # Apply the renaming function to all data frames in the list
    data_list <- lapply(data_list, rename_columns, new_names = new_column_names)
    
    # Specify the desired column order
    new_column_order <- c("deployid","date_time", "longitude", "latitude","t",
                          "x","y")
    
    # Function to reorder columns based on new_column_order
    reorder_columns <- function(df, new_order) {
      df <- df[, new_order]  # Reorder columns by the specified order
      return(df)
    }
    
    # Apply the function to each data frame in the list
    data_list <- lapply(data_list, reorder_columns, 
                        new_order = new_column_order)
    
    collar_data <- rbindlist(data_list)}
  
  if(CSV == TRUE & SHP == FALSE & individual_files == FALSE){
    # Get the list of all CSV files in the directory
    data_file <- read.csv(file, stringsAsFactors = FALSE)
    
    # Define the columns to remove (e.g., "column_name1", "column_name2")
    columns_to_remove <- c("X", "age","dop","collar","collar_deployment_date", 
                           "collar_end_date","duplicate_","fast_step_", 
                           "fast_roundtrip_","nrow")
    
    # Define the new column names
    new_column_names <- c("deployid","x", "y", "date_time","long","lat")
    
    # Function to remove columns from each data frame in the list
    collar_data <- data_file %>% 
      dplyr::select("animal_id",
                    "t_",
                    "longitude",                   
                    "latitude",
                    "x_",                        
                    "y_")  %>% # Select all columns except those in 'cols' 
      rename(., c("deployid" = "animal_id","x" = "x_", "y" = "y_", 
                  "date_time" = "t_", "long" = "longitude", 
                  "lat" = "latitude")) %>% 
      mutate(date_time = ymd_hms(date_time),
             t = as.numeric(lubridate::seconds(date_time))) %>% 
      dplyr::select("deployid","date_time", "long", "lat","t","x","y")
  }
  
  if(SHP == TRUE & CSV == FALSE){
    collar_data <- terra::vect(file)
    
    if("t" %in% names(collar_data))  
    {collar_data$t_ <- as.POSIXct(collar_data$t, 
                                  origin="1970-01-01",
                                  tz="UTC")}
    
    collar_data <- as.data.frame(collar_data)
    collar_data <- collar_data %>% 
      mutate(t_ = parse_date_time(t_, orders = c("ymd_HMS", "ymd")),
             seconds = as.numeric(seconds(t_))) %>% 
      dplyr::select("animal_id","t_","Longitude", "Latitude", "seconds", "UTME",
                    "UTMN") %>%
      rename("deployid" = "animal_id",
             "date_time" = "t_",
             "longitude" = "Longitude", 
             "latitude" = "Latitude",
             "t" = "seconds",
             "x" = "UTME",
             "y" = "UTMN") %>% 
      mutate(longitude = as.numeric(longitude),
             latitude = as.numeric(latitude),
             x = as.numeric(x),
             y = as.numeric(y))
    
  }
  return(collar_data)
}
