# Example code to generate input file for SARS-CoV-2 simulation, using Elm 
# Creek, MN as an example. 

# ECRP_study_area_convex_hull_clipped.shp, MN_metadata.csv, collar_fix_ECRP.RDS,
# EC_County_Parcels.shp, and nlcd_2021_land_cover_metro_area_clipped.tif are not
# included in this repo. This script steps through the creation of a simulation 
# input file, but requires the user to input their own data to generate the 
# resulting file.

# Estimates resource selection for movement (RSF), calculates similarity 
# between observed kernel UDs and grid centroids buffered by mean HR diameter 
# (Origin), estimates exponential decay of use with distance (Decay),
# and prepares parameters for spatial simulations.

# Load Libraries

# Spatial Tools
library(sf)
library(raster)
library(terra)
library(adehabitatHR)

# Tidyverse
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Modeling
library(IndRSA)
library(drc)

### Setup ####

set.seed(567)  # For reproducibility

# Population parameters
pop_size <- NULL # to be assigned at simulation
ave_group_size <- 1.24
sex_ratio <- 28/72

# Load and Prep Location Data

loc_df <- GPS_fix_prep(directory = "~/collar_fix_ECRP.RDS", 
                       start = "2023-09-01", end = "2023-11-30")

annual_loc_df <- GPS_fix_prep(directory = "~/collar_fix_ECRP.RDS", 
                              start = "2023-06-01", end = "2024-05-31")

collar.meta <- collar.meta.prep(filepath = "~/MN_metadata.csv",
                                Site_name = "Elm Creek Park Reserve", 
                                capture_field = "event_type", 
                                capture_status = "First capture", 
                                start_field = "collection_date", 
                                end_field = "deploy_off_date") %>% 
  mutate(animal_id = str_remove(animal_id, "[A-Za-z]$"))


daily_mvmt <- daily_step_length(data = annual_loc_df, CRS = 32615)

median_mvmt <- extract_movement(daily_steps = daily_mvmt, 
                                collar_metadata = collar.meta, 
                                result = "median.mvmt")

mvmt_dist_shape <- extract_movement(daily_steps = daily_mvmt, 
                                    collar_metadata = collar.meta, 
                                    result = "mvmt.dist.shape")

grid <- create_study_area_grid(
  median_mvmt = median_mvmt,
  buffer = 1000,
  common_crs = "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs",
  mask_file = 
    "~/ECRP_study_area_convex_hull_clipped.shp", 
  SHP = TRUE
)

# Create collared animal list for the season. This will be used to randomly 
# assign RSF and origin locations to simulated individuals.

collar_roster <- tibble(animal_id = unique(loc_df$deployid)) %>% 
  left_join(
    collar.meta %>% 
      mutate(animal_id = str_remove(animal_id, "[A-Za-z]$"))
    , by = "animal_id")

###RSF####
#### RSF Part A: Data Prep####

nlcd_rast <- rast(
  "~/nlcd_2021_land_cover_metro_area_clipped.tif")

rcl_mat <- matrix(c(
  11, 1, 21, 2, 22, 3, 23, 4, 24, 5, 31, 5,
  41, 6, 42, 7, 43, 7, 52, 8,
  71, 9, 81, 9, 82, 9,
  90, 10, 95, 10
), ncol = 2, byrow = TRUE)

reclass_rast <- classify(nlcd_rast, rcl = rcl_mat)

grid_RSF <- grid %>% 
  tidyterra::mutate(grid_id = 1:nrow(.))

rast_grid_id <- rasterize(grid_RSF, nlcd_rast, field = "grid_id")
rast_grid_id_cropped <- crop(rast_grid_id, grid_RSF)

# Prep Used vs Unused Points

daily_loc_resampled <- loc_df %>%
  nest(data = -deployid) %>%
  mutate(track = map(data, ~ make_track(tbl = .x, .x = x, .y = y, 
                                        .t = date_time, crs = "+proj=utm 
                                        +zone=15 +datum=WGS84 +units=m 
                                        +no_defs")))

used_unused <- daily_loc_resampled %>% 
  mutate(hr.iso = map(track, ~ hr_mcp(.) |> random_points(n = 10000, 
                                                          presence = .))) %>% 
  mutate(
    covertype = map(hr.iso, ~ extract_covariates(., reclass_rast)),
    grid_id = map(hr.iso, ~ extract_covariates(., rast_grid_id_cropped))
  )

df_used_unused <- map2_dfr(used_unused$covertype, used_unused$deployid, 
                           ~ mutate(as_tibble(.x), deployid = .y)) %>%
  rename(animal_id = deployid)

df_grid_id <- map2_dfr(used_unused$grid_id, used_unused$deployid, 
                       ~ mutate(as_tibble(.x), deployid = .y)) %>%
  rename(animal_id = deployid)

df_used_unused <- bind_cols(df_used_unused, grid_id = df_grid_id$grid_id)
df_used_unused <- left_join(df_used_unused, collar.meta)

# Covariate assignments

assign_cover_type <- function(nlcd_val) {
  case_when(
    nlcd_val == 1 ~ "Water",
    nlcd_val %in% c(2, 3) ~ "open_low_dev",
    nlcd_val %in% c(4, 5) ~ "med_hi_dev",
    nlcd_val %in% c(6, 7, 8) ~ "Forest",
    nlcd_val == 9 ~ "Grass_ag",
    nlcd_val == 10 ~ "Wetlands"
  )
}

df_for_pixel_RSF <- df_used_unused %>% 
  mutate(covertype = assign_cover_type(nlcd_2021_land_cover_metro_area_clipped)
  ) %>% 
  dplyr::select(-nlcd_2021_land_cover_metro_area_clipped) %>% 
  distinct(case_, x_, y_, animal_id, grid_id, covertype, Sex, Age_Class) %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = covertype, values_from = value, values_fill = 
                list(value = 0))

#### RSF Part B: Model Fitting (IndRSA) ####

ls_RSF<-list()
# Fit just one model, as the covariates are binary designations, with 
# open_low_dev as the intercept
ls_RSF[[1]]<-as.formula(case_~Water+Wetlands+Grass_ag+Forest+med_hi_dev)

# Fit individual RSFs
out<-rsf_ind(df_for_pixel_RSF$animal_id, data=df_for_pixel_RSF, 
             form_ls=ls_RSF)

#Extract individual coefficients
pop <- pop_avg(m=1, out, method = "boot")

# Summarize 
ind_selection <- as.data.frame(pop[2]) %>%
  rename("open_low_dev" = "X.Intercept.") %>%
  dplyr::select(-name, -Freq) %>%
  pivot_longer(
    cols = open_low_dev:med_hi_dev,   # Adjust range if needed
    names_to = "covariate",
    values_to = "value"
  ) %>%
  
  # Group by individual animal
  group_by(ID) %>%
  
  # Calculate logit_select_coef relative to open_low_dev category
  mutate(
    open_low_dev_value = value[covariate == "open_low_dev"],
    logit_select_coef = if_else(covariate == "open_low_dev", value, value + 
                                  open_low_dev_value)
  ) %>%
  
  # Scale within animal (using exp-scale normalization)
  mutate(
    logit_select_coef_scaled = exp(logit_select_coef) / 
      max(exp(logit_select_coef), na.rm = TRUE)
  ) %>%
  
  # Select and ungroup final result
  ungroup() %>%
  dplyr::select(ID, covariate, logit_select_coef_scaled) %>%
  mutate(name = "Individual")

#### RSF Part C: Assign selection probabilities across model grid####

# Calculate proportions of covertypes by grid cell
reclass_vect_sample_grid <- terra::extract(reclass_rast, grid_RSF, 
                                           fun = "table", exact = T)

# Convert to proportions and format for use
totals <- reclass_vect_sample_grid %>%
  rowwise() %>%
  mutate(total = sum(c_across(-ID), na.rm = TRUE)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(across(
    -c(ID, total),
    ~ .x / total
  )) %>%
  ungroup() %>%
  dplyr::select(-total) %>%
  rename(., "grid_id" = "ID")

#Convert these proportions back into a SpatVector for prediction
selection_grid <- grid_RSF %>%
  tidyterra::left_join(totals, by = "grid_id") %>%
  tidyterra::select(-"centroid_x", -"centroid_y") %>%
  #tidyterra::filter(!(grid_id %in% water_remove)) %>%
  rowwise() %>%
  mutate(
    Water     = `1`,
    Open_space   = `2`,
    open_low_dev    = sum(c_across(c(`2`, `3`)), na.rm = TRUE),
    med_hi_dev    = sum(c_across(c(`4`, `5`)), na.rm = TRUE),
    Forest    = sum(c_across(c(`6`, `7`, `8`)), na.rm = TRUE),
    Grass_ag  = `9`,
    Wetlands  = `10`
  ) %>%
  ungroup() %>%
  tidyterra::select("grid_id", "Water", "open_low_dev", "med_hi_dev", "Forest", 
                    "Grass_ag", "Wetlands") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

#Calculate aggregated probability on model grid scale
#Rename selection grid
selection_grid_ind_level <- selection_grid

# Convert SpatVector to data.frame
selection_grid_ind_level_df <- as.data.frame(selection_grid_ind_level)

# Convert SpatVector to sf
selection_grid_ind_level_sf <- sf::st_as_sf(selection_grid_ind_level)

# Extract habitat columns (excluding grid_id) to ensure correct column order
habitat_order <- colnames(selection_grid_ind_level_df)[colnames(
  selection_grid_ind_level_df) != "grid_id"]

# Build matrix of grid cell (rows) selection probabilities for each animal 
# (columns)
# Convert to wide matrix: one row per ID, columns = covariates
ind_coef_matrix <- ind_selection %>%
  dplyr::select(ID, covariate, logit_select_coef_scaled) %>%
  pivot_wider(names_from = covariate, 
              values_from = logit_select_coef_scaled) %>%
  arrange(ID)

# Ensure coef_matrix has the same column order (and fill missing with 0)
coef_matrix <- ind_coef_matrix %>%
  dplyr::select(ID, all_of(habitat_order)) %>%
  arrange(ID) %>%
  rowwise() %>%
  mutate(across(
    all_of(habitat_order),
    ~ if_else(is.na(.x), open_low_dev, .x)
  )) %>%
  ungroup() %>%
  dplyr::select(-ID)

# Convert grid to matrix for multiplication (remove grid_id temporarily)
grid_matrix <- selection_grid_ind_level_df %>%
  dplyr::select(all_of(habitat_order)) %>%
  as.matrix()

# Matrix multiplication
selection_matrix <- grid_matrix %*% t(coef_matrix)

#Create final dataframe for RSF in simulation, 
# setting grid_ids with 100% water to 0
selection_grid_with_individuals <- selection_grid_ind_level_df %>%
  bind_cols(as.data.frame(selection_matrix) %>% 
              setNames(ind_coef_matrix$ID)) %>% 
  dplyr::select(!all_of(habitat_order)) %>% 
  mutate(across(all_of(ind_coef_matrix$ID),
                ~ if_else(grid_id %in% selection_grid_ind_level_df[
                  which(selection_grid_ind_level_df$Water == 1),"grid_id"], 0,
                  .x)))

# Outputs for this section used in simulations 
# include: selection_grid_with_individuals

### Origin ####
#### Origin Part A: Data Prep ####

# Arrange data for kernel UD analysis
loc_by_id <- loc_df %>% 
  dplyr::select(animal_id = deployid, date_time, x, y) %>% 
  merge(., collar.meta, by = c("deploy_id" = "animal_id")) %>% 
  st_as_sf(., coords = 3:4, crs = "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs")

# Subset female individuals
female_loc <- loc_by_id %>% 
  filter(Sex == "F")

# Subset male individuals
male_loc <- loc_by_id %>% 
  filter(Sex == "M")


#### Origin Part B: Kernel UD calculation and habitat composition ####


# Calculate UDs for female individuals 
female_deer_ud <- adehabitatHR::kernelUD(as(female_loc, "Spatial")[1],
                                         grid = 500, extent = 1)

# Calculate UDs for male individuals
male_deer_ud <- adehabitatHR::kernelUD(as(male_loc, "Spatial")[1],
                                       grid = 500, extent = 1)

# Estimate size of UD to establish sampling radius 
# (presuming circular home range)
mean_female_radius <- sqrt(mean(as.numeric(adehabitatHR::kernel.area(
  female_deer_ud, percent = 95, unin = "m", unout = "m")))/pi)

mean_male_radius <- sqrt(mean(as.numeric(adehabitatHR::kernel.area(
  male_deer_ud, percent = 95, unin = "m", unout = "m")))/pi)

# Function to extract proportions
extract_ud_habitat_props <- function(ud, habitat_raster, animal_id) {
  # Convert UD to raster (if not already)
  ud_raster <- raster(ud)
  
  # Resample UD to match habitat resolution and extent (if needed)
  ud_resampled <- resample(ud_raster, habitat_raster, method = "bilinear")
  
  # Mask habitat raster by UD (optional but speeds up computation)
  masked_habitat <- raster::mask(habitat_raster, ud_resampled)
  
  # Extract values
  habitat_vals <- getValues(masked_habitat)
  ud_vals <- getValues(ud_resampled)
  
  # Remove missing data
  valid_idx <- which(!is.na(habitat_vals) & !is.na(ud_vals))
  habitat_vals <- habitat_vals[valid_idx]
  ud_vals <- ud_vals[valid_idx]
  
  # Compute weighted sum of habitat proportions
  prop_df <- tibble(habitat = habitat_vals, weight = ud_vals) %>%
    group_by(habitat) %>%
    summarise(prop = sum(weight, na.rm = TRUE)) %>%
    mutate(prop = prop / sum(prop)) %>%  # normalize to sum to 1
    mutate(ID = animal_id)
  
  return(prop_df)
}

# Save the reclass raster (terra package) to raster (raster package) for 
# function to work.
reclass_raster <- raster(reclass_rast)

# Apply over all females
female_habitat_props_list <- imap(female_deer_ud, 
                                  ~ extract_ud_habitat_props(.x, reclass_raster,
                                                             .y))

# Combine to one dataframe
female_habitat_props_df <- bind_rows(female_habitat_props_list) %>% 
  mutate(habitat = factor(case_when(habitat == 1 ~ "Water",
                                    habitat %in% c(2,3) ~ "open_low_dev",
                                    habitat %in% c(4,5) ~ "med_hi_dev",
                                    habitat %in% c(6,7,8) ~ "Forest",
                                    habitat == 9 ~ "Grass_ag",
                                    habitat == 10 ~ "Wetlands"), 
                          levels = habitat_order)) %>% 
  group_by(ID, habitat) %>% 
  reframe(prop = sum(prop)) %>% 
  ungroup()

# Prep habitat composition for use
female_habitat_props_df <- female_habitat_props_df %>% 
  pivot_wider(names_from = habitat, values_from = prop, 
              values_fill = 0) %>% 
  dplyr::select(c(ID, all_of(habitat_order)))

# Apply over all males
male_habitat_props_list <- imap(male_deer_ud, 
                                ~ extract_ud_habitat_props(.x, reclass_raster, 
                                                           .y))

# Combine to one dataframe
male_habitat_props_df <- bind_rows(male_habitat_props_list) %>% 
  mutate(habitat = factor(case_when(habitat == 1 ~ "Water",
                                    habitat %in% c(2,3) ~ "open_low_dev",
                                    habitat %in% c(4,5) ~ "med_hi_dev",
                                    habitat %in% c(6,7,8) ~ "Forest",
                                    habitat == 9 ~ "Grass_ag",
                                    habitat == 10 ~ "Wetlands"), 
                          levels = habitat_order)) %>% 
  group_by(ID, habitat) %>% 
  reframe(prop = sum(prop)) %>% 
  ungroup()

# Prep habitat composition for use
male_habitat_props_df <- male_habitat_props_df %>% 
  pivot_wider(names_from = habitat, values_from = prop, 
              values_fill = 0) %>% 
  dplyr::select(c(ID, all_of(habitat_order)))


#### Origin Part C: Compare habitat composition around each grid cell ####

# Establish grid cell centroids
grid_points <- centroids(grid)

# Buffer centroids by HR radius, by sex
grid_point_buffered_female <- terra::buffer(grid_points, width = mean_female_radius)
grid_point_buffered_male <- terra::buffer(grid_points, width = mean_male_radius)

# Intersect buffers to study area extent
grid_union <- terra::aggregate(grid)

grid_point_buffered_female_clipped <- terra::intersect(grid_point_buffered_female, grid_union)
grid_point_buffered_male_clipped <- terra::intersect(grid_point_buffered_male, grid_union)

# Extract habitat values from raster within each buffer
habitat_extract_female <- terra::extract(reclass_rast, grid_point_buffered_female_clipped, fun = NULL, cells = FALSE)
habitat_extract_male <- terra::extract(reclass_rast, grid_point_buffered_male_clipped, fun = NULL, cells = FALSE)

# Combine with grid IDs
habitat_extract_female$grid_id <- grid_point_buffered_female_clipped$grid_id[habitat_extract_female$ID]
habitat_extract_male$grid_id <- grid_point_buffered_male_clipped$grid_id[habitat_extract_male$ID]

# Extract habitat composition for each grid cell centroid buffered by 
# female radius
habitat_props_grid_female <- habitat_extract_female %>%
  mutate(
    habitat = factor(
      case_when(
        nlcd_2021_land_cover_metro_area_clipped == 1 ~ "Water",
        nlcd_2021_land_cover_metro_area_clipped %in% c(2, 3) ~ "open_low_dev",
        nlcd_2021_land_cover_metro_area_clipped %in% c(4, 5) ~ "med_hi_dev",
        nlcd_2021_land_cover_metro_area_clipped %in% c(6, 7, 8) ~ "Forest",
        nlcd_2021_land_cover_metro_area_clipped == 9 ~ "Grass_ag",
        nlcd_2021_land_cover_metro_area_clipped == 10 ~ "Wetlands"
      ),
      levels = habitat_order
    )
  ) %>%
  filter(!is.na(habitat)) %>%
  group_by(grid_id, habitat) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(grid_id) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = grid_id,
    names_from = habitat,
    values_from = prop,
    values_fill = 0
  ) %>%
  dplyr::select(grid_id, all_of(habitat_order))

# Extract habitat composition for each grid cell centroid buffered by 
# male radius
habitat_props_grid_male <- habitat_extract_male %>%
  mutate(
    habitat = factor(
      case_when(
        nlcd_2021_land_cover_metro_area_clipped == 1 ~ "Water",
        nlcd_2021_land_cover_metro_area_clipped %in% c(2, 3) ~ "open_low_dev",
        nlcd_2021_land_cover_metro_area_clipped %in% c(4, 5) ~ "med_hi_dev",
        nlcd_2021_land_cover_metro_area_clipped %in% c(6, 7, 8) ~ "Forest",
        nlcd_2021_land_cover_metro_area_clipped == 9 ~ "Grass_ag",
        nlcd_2021_land_cover_metro_area_clipped == 10 ~ "Wetlands"
      ),
      levels = habitat_order
    )
  ) %>%
  filter(!is.na(habitat)) %>%
  group_by(grid_id, habitat) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(grid_id) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = grid_id,
    names_from = habitat,
    values_from = prop,
    values_fill = 0
  ) %>%
  dplyr::select(grid_id, all_of(habitat_order))


#### Origin Part D: Calculate similarity between buffered grid cell and ud####

# Female comparison
# Match habitat column order
common_habitats_female <- intersect(names(habitat_props_grid_female), names(female_habitat_props_df))

# Ensure proper habitat column order and arrange by grid_id for grid matrix
female_grid_matrix <- habitat_props_grid_female %>%
  dplyr::select(grid_id, all_of(common_habitats_female)) %>%
  arrange(grid_id)

# Ensure proper habitat column order and arrange by individual for 
# individual matrix
female_ind_matrix <- female_habitat_props_df %>%
  dplyr::select(ID, all_of(common_habitats_female)) %>%
  arrange(ID)

# Convert ind_matrix to a named list of numeric vectors (one per individual)
female_ind_vectors <- female_ind_matrix %>%
  column_to_rownames("ID") %>%
  split(factor(rownames(.))) %>%  # Split by row
  map(~ as.numeric(.x))    # Convert each to numeric vector

# Compute similarity correctly, as inverse of sum of absolute value difference 
# between grid habitat availability and UD habitat availability. I add 1 to the 
# denominator so that the result is bounded between 0 and 1.
female_similarity_matrix <- female_grid_matrix %>%
  rowwise() %>%
  mutate(
    similarities = list(
      map_dbl(female_ind_vectors, ~ 1 / 
                (sum(abs(c_across(all_of(common_habitats_female)) - .x))+1))
    )
  ) %>%
  ungroup()

#Find global min and max
all_values <- unlist(female_similarity_matrix$similarities)
min_val <- min(all_values)
max_val <- max(all_values)

#Scale selection value for each animal
female_similarity_matrix <- female_similarity_matrix %>%
  mutate(similarities_scaled = map(similarities, 
                                   ~ (.-min_val)/(max_val-min_val)))

# Convert similarity list column into wide format
female_similarity_df <- female_similarity_matrix %>%
  dplyr::select(grid_id, similarities_scaled) %>%
  unnest_wider(similarities_scaled, names_sep = "_") %>%
  rename_with(~ names(female_ind_vectors), 
              starts_with("similarities_scaled_")) %>% 
  mutate(across(all_of(female_ind_matrix$ID),
                ~ if_else(grid_id %in% selection_grid_ind_level_df[which(selection_grid_ind_level_df$Water == 1),"grid_id"], 0, .x)))

# Male comparison
# Match habitat column order
common_habitats_male <- intersect(names(habitat_props_grid_male), names(male_habitat_props_df))

# Ensure proper habitat column order and arrange by grid_id for grid matrix
male_grid_matrix <- habitat_props_grid_male %>%
  dplyr::select(grid_id, all_of(common_habitats_male)) %>%
  arrange(grid_id)

# Ensure proper habitat column order and arrange by individual for 
# individual matrix
male_ind_matrix <- male_habitat_props_df %>%
  dplyr::select(ID, all_of(common_habitats_male)) %>%
  arrange(ID)

# Convert ind_matrix to a named list of numeric vectors (one per individual)
male_ind_vectors <- male_ind_matrix %>%
  column_to_rownames("ID") %>%
  split(factor(rownames(.))) %>%  # Split by row
  map(~ as.numeric(.x))    # Convert each to numeric vector

# Compute similarity correctly, as inverse of sum of absolute value difference 
# between grid habitat availability and UD habitat availability. I add 1 to the 
# denominator so that the result is bounded between 0 and 1.
male_similarity_matrix <- male_grid_matrix %>%
  rowwise() %>%
  mutate(
    similarities = list(
      map_dbl(male_ind_vectors, ~ 1 / 
                (sum(abs(c_across(all_of(common_habitats_male)) - .x))+1))
    )
  ) %>%
  ungroup()
#Find global min and max
all_values <- unlist(male_similarity_matrix$similarities)
min_val <- min(all_values)
max_val <- max(all_values)

#Scale selection value for each animal
male_similarity_matrix <- male_similarity_matrix %>%
  mutate(similarities_scaled = map(similarities, 
                                   ~ (.-min_val)/(max_val-min_val)))

# Convert similarity list column into wide format
male_similarity_df <- male_similarity_matrix %>%
  dplyr::select(grid_id, similarities_scaled) %>%
  unnest_wider(similarities_scaled, names_sep = "_") %>%
  rename_with(~ names(male_ind_vectors), 
              starts_with("similarities_scaled_")) %>% 
  mutate(across(all_of(male_ind_matrix$ID),
                ~ if_else(grid_id %in% selection_grid_ind_level_df[which(selection_grid_ind_level_df$Water == 1),"grid_id"], 0, .x)))

similarity_df <- left_join(female_similarity_df, male_similarity_df)

# Outputs for this section used in simulations include: female_similarity_df 
# and male_similarity_df

### Centroid fidelity decay####
#### Decay Part A: Identify location of greatest utilization (centroid)####

# Create UD for all collared animals, add the id column and you get unique 
# homeranges
deer_ud <- adehabitatHR::kernelUD(as(loc_by_id, "Spatial")[1],
                                  grid = 500, extent = 1)

# Extract home ranges for each animal
hr_vol <- adehabitatHR::getvolumeUD(deer_ud, standardize = T)

# Convert UD to raster
ud_rasters <- lapply(1:length(hr_vol), function(i) {
  raster(as(hr_vol[[i]], "SpatialPixelsDataFrame"))
})
names(ud_rasters) <- names(hr_vol)  # Assign animal IDs as names

# Extract location of raster pixel with greatest use (minimum value)
centroids <- lapply(ud_rasters, function(r) {
  coords <- xyFromCell(r, 1:ncell(r))  # Get coordinates of all pixels
  values <- getValues(r)  # Extract UD values
  
  centroid <- coords[which(values==min(values)),]
  return(centroid)
})


# Convert to a data frame for easy reference
centroid_df <- as_tibble(do.call(rbind, centroids)) %>% 
  #rename(., "x" = V1, "y" = V2) %>% 
  mutate(animal_id = names(ud_rasters))

#Convert centroid locations into SpatVect (terra package)
centroid_vect <- vect(centroid_df, geom = c("x", "y"), crs = "+proj=utm +zone=15
                      +datum=WGS84 +units=m +no_defs")

# Also save as an sf object
centroid_sf <- st_as_sf(centroid_vect)


#### Decay Part B: Calculate distances from centroid for observed and available#### 

# Check CRS consistency
if (crs(loc_by_id) != crs(centroid_vect)) {
  stop("CRS of loc_by_id and centroid_vect do not match. Reproject if needed.")
} else {
  message("CRS of loc_by_id and centroid_vect match.")
}

# Convert SpatVectors to data frames and extract coordinates
centroid_df <- as.data.frame(centroid_vect, geom = "XY") #

# Bring in samples used for RSF section; rename the coordinate columns 
# for clarity
fixes <- df_for_pixel_RSF %>%
  rename(x_loc = x_, y_loc = y_)

# Format centroid dataframe with correct order
centroid_df <- centroid_df %>%
  rename(x_centroid = x, y_centroid = y) %>% 
  dplyr::select(animal_id, x_centroid, y_centroid)

# Select fixes with the relevant columns for merging
fixes <- fixes %>%
  dplyr::select(case_, animal_id, Sex, Age_Class, x_loc, y_loc)

# Merge the location data with the centroid data based on animal_id
merged_data <- left_join(fixes, centroid_df, by = "animal_id")

# Calculate the distance using the coordinates for locations and centroid
merged_data <- merged_data %>%
  mutate(distance_to_centroid = sqrt((x_loc - x_centroid)^2 + 
                                       (y_loc - y_centroid)^2))

#Create complete dataframe of distance of locations to centroids
df_for_distance_decay <- merged_data %>% 
  dplyr::select(case_, animal_id, Sex, Age_Class, distance_to_centroid) %>% 
  mutate(animal_id = as.factor(animal_id),
         Sex = as.factor(Sex),
         Age_Class = as.factor(Age_Class),
         distance_to_centroid = distance_to_centroid/1000)

# Reduce the number of unused points from 10,000 to a number equal to 
# the used points
df_equal <- df_for_distance_decay %>%
  group_by(animal_id) %>%
  nest() %>%
  mutate(
    true_rows = map(data, ~ filter(.x, case_ == TRUE)),
    false_rows = map(data, ~ filter(.x, case_ == FALSE)),
    n_true = map_int(true_rows, nrow),
    min_count = min(n_true),
    false_sample = map2(false_rows, n_true, ~ {
      n_sample <- min(nrow(.x), .y)
      sample_n(.x, n_sample)
    })
  ) %>%
  dplyr::select(animal_id, true_rows, false_sample) %>%
  unnest(c(true_rows, false_sample), names_sep = "_") %>%
  ungroup()

# Create dataframe for observed locations (T) and random locations (F)
df_equal_T <- df_equal[,c("animal_id", "true_rows_Sex", "true_rows_Age_Class", "true_rows_distance_to_centroid")] %>% 
  rename(., "Sex" = "true_rows_Sex", "Age_Class" = "true_rows_Age_Class", "distance_to_centroid" = "true_rows_distance_to_centroid") %>% 
  mutate(case_ = TRUE)

df_equal_F <- df_equal[,c("animal_id", "false_sample_Sex", 
                          "false_sample_Age_Class",
                          "false_sample_distance_to_centroid")] %>% 
  rename(., "Sex" = "false_sample_Sex", "Age_Class" = "false_sample_Age_Class", "distance_to_centroid" = "false_sample_distance_to_centroid") %>% 
  mutate(case_ = FALSE)

# And combine
df_equal <- rbind(df_equal_T, df_equal_F)

#### Decay Part C: Fit exponential decay models as function of distance, by sex####

# Isolate female individuals
female_df_nls <- df_equal %>% 
  filter(., Sex =="F")

# Fit exponential decay function for females
dist_mod_drs_f_exp <- drm(case_ ~ distance_to_centroid, type = "binomial", 
                          fct = EXD.2(fixed = c(1,NA)),
                          data = female_df_nls)

# Predict use as a function of distance
pref_distance_F <- data.frame(distance_to_centroid = seq(0, max(female_df_nls$distance_to_centroid), length.out = 50)) %>% 
  mutate(pred_use_exp = predict(dist_mod_drs_f_exp, .),
         Sex = "F")

# Isolate male individuals
male_df_nls <- df_equal %>% 
  filter(., Sex =="M") 

# Fit exponential decay function for males
dist_mod_drs_m_exp <- drm(case_ ~ distance_to_centroid, type = "binomial", 
                          fct = EXD.2(fixed = c(1,NA)),
                          data = male_df_nls)

#Predict using female max distance
pref_distance_M <- data.frame(distance_to_centroid = seq(0, max(female_df_nls$distance_to_centroid), length.out = 50)) %>% 
  mutate(pred_use_exp = predict(dist_mod_drs_m_exp, .),
         Sex = "M")

pref_distance <- rbind(pref_distance_F, pref_distance_M)

# Plotting code for distance to centroid
exp_plot <- ggplot(pref_distance, aes(x = distance_to_centroid, 
                                      y = pred_use_exp, color = Sex))+
  geom_line()+
  theme_classic() +
  labs(tag = "2-parameter exponential decay")+
  scale_x_continuous(name = "Distance from centroid (km)") +
  scale_y_continuous(name = "Proportion of cells used by individuals")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

# Package decay with distance results for modeling, storing all predictions in 
# one dataframe for reference during simulations
centroid_distance_decay <- list(female_centroid_distance_decay = 
                                  dist_mod_drs_f_exp, 
                                male_centroid_distance_decay = 
                                  dist_mod_drs_m_exp)

# Check that individuals are consistent across objects, and remove any that are

# Compare 3rd order RSF with roster (and vice versa). This code will identify if
# an individual is listed in one object but not the other.
colnames(selection_grid_with_individuals[-1])[which(!(colnames(selection_grid_with_individuals[,-1]) %in% collar_roster$animal_id))]
collar_roster$animal_id[which(!(collar_roster$animal_id %in% colnames(selection_grid_with_individuals[,-1])))]

# Compare 2nd order RSF with roster (and vice versa): This code will identify if an individual is listed in one object but not the other.
colnames(similarity_df[-1])[which(!(colnames(similarity_df[,-1]) %in% collar_roster$animal_id))]
collar_roster$animal_id[which(!(collar_roster$animal_id %in% colnames(similarity_df[,-1])))]

# Compare 2nd order to 3rd order 
colnames(similarity_df[-1])[which(!(colnames(similarity_df[,-1]) %in% colnames(selection_grid_with_individuals[,-1])))]
colnames(selection_grid_with_individuals[-1])[which(!(colnames(selection_grid_with_individuals[,-1]) %in% colnames(similarity_df[,-1])))]

### Export!####


sim_input <- list(pop_size = pop_size, ave_group_size = ave_group_size, sex_ratio = sex_ratio, collar_roster = collar_roster, selection_probs = selection_grid_with_individuals, origin_probs = similarity_df, centroid_decay = centroid_distance_decay)

saveRDS(sim_input, "~/ECRP_sim_input_fall.RDS")
