# Create contact rate - distance relationship, using Elm Creek, MN as example. 
# This same script was applied to the other 8 study areas. Collar location data 
# is not included in this repo, so this script is for demonstration purposes.

####
#ECRP
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/MN_metadata.csv",
                                Site_name = "Elm Creek Park Reserve", 
                                capture_field = "event_type", 
                                capture_status = "First capture",
                                start_field = "collection_date", 
                                end_field = "deploy_off_date")

ECRP_annual_loc_df <- GPS_fix_prep(directory = "~/collar_fix_ECRP.RDS", 
                                   start = "2023-06-01", end = "2024-05-31")

ECRP_annual_dist_daily.df <- daily_centroid_dist(animal_data = 
                                                   ECRP_annual_loc_df, 
                                                 collar_meta = collar.meta, 
                                                 centroid_crs = '+proj=utm 
                                                 +zone=15 ellps=WGS84', 
                                                 grid_size = 30, 
                                                 extent_margin = 500, 
                                                 min_points_per_day = 5, 
                                                 num_cores = 64)

ECRP.contact.est.annual <- calc_contact_rates_continuous_gamma(
                              collar_locations = ECRP_annual_loc_df, 
                              contact_crs = '+proj=utm +zone=15 ellps=WGS84',
                              fix_animal_field = "deployid", 
                              time_field = "date_time",
                              distance_tol = 25, 
                              time_tol = 15*60, 
                              collar_metadata = collar.meta, 
                              collar_metadata_animal_field =  "animal_id", 
                              centroid_distances = ECRP_annual_dist_daily.df, 
                              start.date = "2023-06-01",
                              end.date = "2024-05-31", num_cores = 64)

#saveRDS(ECRP.contact.est.annual, file = "~/contact_rates_annual_EC.RDS")
