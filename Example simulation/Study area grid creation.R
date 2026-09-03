## Code to delineate a study area, using Elm Creek, MN as an example. 
# ECRP_study_area_convex_hull_clipped.shp, MN_metadata.csv, collar_fix_ECRP.RDS,
# and EC_County_Parcels.shp are not included in this repo.

#Bring in buffered study area polygon
ECRP_study_area <-  vect("~/ECRP_study_area_convex_hull_clipped.shp", 
                         crs = "+proj=utm +zone=15 +datum=WGS84 
                         +units=m +no_defs +type=crs")

#Read in collar metadata
ECRP.collar.meta <- collar.meta.prep(filepath = "~/MN_metadata.csv",
                                     Site_name = "Elm Creek Park Reserve", 
                                     capture_field = "event_type", 
                                     capture_status = "First capture", 
                                     start_field = "collection_date", 
                                     end_field = "deploy_off_date")

#Bring in annual locations for grid development
ECRP_annual_loc_df <- GPS_fix_prep(directory = "~/collar_fix_ECRP.RDS", 
                                   start = "2023-06-01", end = "2024-05-31")

#Extract daily movements
ECRP_daily_mvmt <- daily_step_length(data = ECRP_annual_loc_df, 
                                     CRS = "+proj=utm +zone=15 +datum=WGS84 
                                     +units=m +no_defs +type=crs")

#Extract median daily distance moved by sex/age class
ECRP_median_mvmt <- extract_movement(daily_steps = ECRP_daily_mvmt, 
                                     collar_metadata = ECRP.collar.meta, 
                                     result = "median.mvmt")

#Create grid for spatial modeling as a spatVector
ECRP.grid <- create_study_area_grid(median_mvmt = ECRP_median_mvmt, 
                                    buffer = 1000, 
                                    common_crs = "+proj=utm +zone=15 
                                    +datum=WGS84 +units=m +no_defs +type=crs", 
                                    mask_file = 
                                    "~/ECRP_study_area_convex_hull_clipped.shp", 
                                    SHP = TRUE)

ECRP.grid <- add_centroid_count_to_grid(
  shapefile_path = "~/EC_County_Parcels.shp", 
  name = "Residence_count",
  grid = ECRP.grid, 
  crs = "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +type=crs")

