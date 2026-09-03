#Create and compare contact rate curves

####
#ECRP
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/MN_Sampling_MasterSheet_use.csv",Site_name = "Elm Creek Park Reserve", 
                                capture_field = "event_type", capture_status = "First capture", start_field = "collection_date", end_field = "deploy_off_date")

ECRP_annual_loc_df <- GPS_fix_prep(directory = "~/ER_Data_Files/collar_data/ECPR/collar_fix_ECRP.RDS", start = "2023-06-01", end = "2024-05-31")

ECRP_annual_dist_daily.df <- daily_centroid_dist(animal_data = ECRP_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=15 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

ECRP.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = ECRP_annual_loc_df, 
                                                             contact_crs = '+proj=utm +zone=15 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                             collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = ECRP_annual_dist_daily.df, start.date = "2023-06-01",
                                                             end.date = "2024-05-31", num_cores = 64)

saveRDS(ECRP.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/ECRP/contact_rates_annual_EC.RDS")

####
#CRP
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/MN_Sampling_MasterSheet_use.csv",Site_name = "Carver Park Reserve", 
                                capture_field = "event_type", capture_status = "First capture", start_field = "collection_date", end_field = "deploy_off_date")

carver_loc_df <- read_raw_GPS_fixes(file = 
                                      "~/ER_Data_Files/collar_data/Carver/Carver_collars.shp", 
                                    SHP = TRUE)

CRP_annual_loc_df <- GPS_fix_prep(vector_object = carver_loc_df, start = "2024-03-01",
                                  end = "2025-02-28")

CRP_annual_dist_daily.df <- daily_centroid_dist(animal_data = CRP_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=15 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

CRP.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = CRP_annual_loc_df, 
                                                               contact_crs = '+proj=utm +zone=15 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                               collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = CRP_annual_dist_daily.df, start.date = "2024-03-01",
                                                               end.date = "2025-02-28", num_cores = 64)

saveRDS(CRP.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/Carver/contact_rates_annual_CP.RDS")

####
#Shakopee
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/MN_Sampling_MasterSheet_use.csv",Site_name = "Shakopee Mdewakanton Sioux Community", 
                                capture_field = "event_type", capture_status = "First capture", start_field = "collection_date", end_field = "deploy_off_date")

#Prep data points for the season
shakopee_loc_df <- read_raw_GPS_fixes(file = 
                                      "~/ER_Data_Files/collar_data/Shakopee/shakopee_collars.shp", 
                                    SHP = TRUE)

shakopee_annual_loc_df <- GPS_fix_prep(vector_object = shakopee_loc_df, start = "2024-03-01",
                                  end = "2025-02-28")

shakopee_annual_dist_daily.df <- daily_centroid_dist(animal_data = shakopee_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=15 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

shakopee.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = shakopee_annual_loc_df, 
                                                              contact_crs = '+proj=utm +zone=15 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                              collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = shakopee_annual_dist_daily.df, start.date = "2024-03-01",
                                                              end.date = "2025-02-28", num_cores = 64)

saveRDS(shakopee.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/Shakopee/contact_rates_annual_shakopee.RDS")

####
#TON - Touch of Nature, IL
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/IL_deer_meta_v04.30.2025.csv",Site_name = "Touch of Nature", 
                                capture_field = "event_type", capture_status = "New", start_field = "Capture_date", end_field = "Deploy_off")

TON_loc_df <- read_raw_GPS_fixes(directory = "~/ER_Data_Files/collar_data/TON", SHP =FALSE, CSV = TRUE, individual_files = TRUE)

TON_annual_loc_df <- GPS_fix_prep(vector_object = TON_loc_df, start = "2023-03-01",
                                  end = "2024-02-29")

TON_annual_dist_daily.df <- daily_centroid_dist(animal_data = TON_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=16 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

TON.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = TON_annual_loc_df, 
                                                              contact_crs = '+proj=utm +zone=16 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                              collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = TON_annual_dist_daily.df, start.date = "2023-03-01",
                                                              end.date = "2024-02-29", num_cores = 64)

saveRDS(TON.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/TON/contact_rates_annual_TON.RDS")

####
#Shelbyville - Shelbyville, IL (SIU)
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/IL_deer_meta_v04.30.2025.csv",Site_name = "Shelbyville", 
                                capture_field = "event_type", capture_status = "New", start_field = "Capture_date", end_field = "Deploy_off")

Shelby_loc_df <- read_raw_GPS_fixes(directory = "~/ER_Data_Files/collar_data/Shelbyville", SHP =FALSE, CSV = TRUE, individual_files = TRUE)

Shelby_annual_loc_df <- GPS_fix_prep(vector_object = Shelby_loc_df, start = "2024-03-01",
                                  end = "2025-02-28")

Shelby_annual_dist_daily.df <- daily_centroid_dist(animal_data = Shelby_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=16 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

Shelby.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = Shelby_annual_loc_df, 
                                                              contact_crs = '+proj=utm +zone=16 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                              collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = Shelby_annual_dist_daily.df, start.date = "2024-03-01",
                                                              end.date = "2025-02-28", num_cores = 64)

saveRDS(Shelby.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/Shelbyville/contact_rates_annual_Shelby.RDS")


####
#Ames - Ames Plantation, TN (UTK)
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/utk_collar_metadata_updated.csv",Site_name = "Ames", 
                                capture_field = "collar", capture_status = "Yes", start_field = "Date", end_field = "Removed_from_study_date")

ames_loc_df <- read_raw_GPS_fixes(file = "~/ER_Data_Files/collar_data/Ames/Ames_collars.shp", SHP =TRUE, CSV = FALSE, individual_files = FALSE)

ames_annual_loc_df <- GPS_fix_prep(vector_object = ames_loc_df, start = "2024-03-01",
                                     end = "2025-02-28")

ames_annual_dist_daily.df <- daily_centroid_dist(animal_data = ames_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=16 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

ames.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = ames_annual_loc_df, 
                                                                 contact_crs = '+proj=utm +zone=16 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                                 collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = ames_annual_dist_daily.df, start.date = "2024-03-01",
                                                                 end.date = "2025-02-28", num_cores = 64)

saveRDS(ames.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/Ames/contact_rates_annual_ames.RDS")

####
#Lone Oaks - Lone Oaks, TN (UTK)
#Prep metadata (Ames is the site for all collars from TN)
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/utk_collar_metadata_updated.csv",Site_name = "Ames", 
                                capture_field = "collar", capture_status = "Yes", start_field = "Date", end_field = "Removed_from_study_date")

lone_oaks_loc_df <- read_raw_GPS_fixes(file = "~/ER_Data_Files/collar_data/Lone_oaks/lone_oaks_collars.shp", SHP =TRUE, CSV = FALSE, individual_files = FALSE)

lone_oaks_annual_loc_df <- GPS_fix_prep(vector_object = lone_oaks_loc_df, start = "2024-09-01",
                                   end = "2025-08-31")

lone_oaks_annual_dist_daily.df <- daily_centroid_dist(animal_data = lone_oaks_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=16 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

lone_oaks.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = lone_oaks_annual_loc_df, 
                                                               contact_crs = '+proj=utm +zone=16 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                               collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = lone_oaks_annual_dist_daily.df, start.date = "2024-09-01",
                                                               end.date = "2025-08-31", num_cores = 64)

saveRDS(lone_oaks.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/Lone_oaks/contact_rates_annual_lone_oaks.RDS")

####
#Staten Island, NY (Columbia)
#Prep metadata
collar.meta <- collar.meta.prep(filepath = "~/ER_Data_Files/collar_data/Odocoileus virginianus Staten Island, NY-reference-data_ER.csv",Site_name = "SI", 
                                capture_field = "site_name", capture_status = "SI", start_field = "deploy_on_date", end_field = "animal_mortality_date")

SI_loc_df <- read_raw_GPS_fixes(file = 
                                    "~/ER_Data_Files/collar_data/SI/SI_collars_combo_clipped.shp", 
                                  SHP = TRUE)

SI_annual_loc_df <- GPS_fix_prep(vector_object = SI_loc_df, start = "2024-03-01",
                                        end = "2025-02-28")

SI_annual_dist_daily.df <- daily_centroid_dist(animal_data = SI_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=18 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

SI.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = SI_annual_loc_df, 
                                                                    contact_crs = '+proj=utm +zone=18 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                                    collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = SI_annual_dist_daily.df, start.date = "2024-03-01",
                                                                    end.date = "2025-02-28", num_cores = 64)

saveRDS(SI.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/SI/contact_rates_annual_SI.RDS")

####
# Treasure Lake, PA
collar.meta <- read.csv(file = "ER_Data_Files/collar_data/TL_collar_metadata.csv") %>% 
  rename("Sex" = "sex",
         "start_deployment" = "collar_deployment_date",
         "end_deployment" = "collar_end_date") %>% 
  mutate(Age_Class = "NA") 


#Prep data points for the season
TL_loc_df <- read_raw_GPS_fixes(file = 
                                    "~/ER_Data_Files/collar_data/TL/PA_TL_collars.shp", 
                                  SHP = TRUE)

TL_annual_loc_df <- GPS_fix_prep(vector_object = TL_loc_df, start = "2024-01-26",
                                   end = "2024-09-02")

TL_annual_dist_daily.df <- daily_centroid_dist(animal_data = TL_annual_loc_df, collar_meta = collar.meta, centroid_crs = '+proj=utm +zone=17 ellps=WGS84', grid_size = 30, extent_margin = 500, min_points_per_day = 5, num_cores = 64)

TL.contact.est.annual <- calc_contact_rates_continuous_gamma(collar_locations = TL_annual_loc_df, 
                                                               contact_crs = '+proj=utm +zone=17 ellps=WGS84',fix_animal_field = "deployid", time_field = "date_time",distance_tol = 25, time_tol = 15*60, 
                                                               collar_metadata = collar.meta, collar_metadata_animal_field =  "animal_id", centroid_distances = TL_annual_dist_daily.df, start.date = "2024-01-26",
                                                               end.date = "2024-09-02", num_cores = 64)

saveRDS(TL.contact.est.annual, file = "~/Results/Contact_rate_selection/contact_rate_functions/TL/contact_rates_annual_TL.RDS")
