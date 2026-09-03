## Elm Creek Reserve Park (ECRP) inputs ####

ECRP_fall_sim_input <- readRDS("Example simulation/Inputs/ECRP_sim_input_fall.RDS")

# See Study area grid creation.R for illustration of building ECRP.grid, 
# ECRP_median_mvmt, and ECRP_mvmt_dist_shape
ECRP.grid <- readRDS("Example simulation/Inputs/ECRP_grid.RDS")

ECRP_median_mvmt <- readRDS("Example simulation/Inputs/ECRP_median_mvmt.RDS")

ECRP_mvmt_dist_shape <- readRDS("Example simulation/Inputs/ECRP_mvmt_shape_params.RDS")

#Create grid dataframe for simulations
ECRP_centroids_fall_df <- ECRP.grid %>% 
  tidyterra::select(centroid_x, centroid_y, grid_id) %>% 
  as.data.frame() %>% 
  rename(c("X" = "centroid_x", "Y" = "centroid_y")) %>% 
  left_join(ECRP_fall_sim_input$selection_probs, by = "grid_id")

# Create a dataframe with pair-wise pixel distances, for use in subsetting all
# potential destination pixels (pixel 2) from an origin pixel (pixel 1), within 
# a reasonable distance for adult males to travel. This cuts down computaion 
# time.
ECRP_fall_dist_btw_pixels <- pixel_options(grid = ECRP.grid, 
                                           movement_params = ECRP_mvmt_dist_shape, 
                                           age_sex = "Adult male", threshold = 0.99, 
                                           centroid_decay_model = ECRP_fall_sim_input$centroid_decay)

## Common inputs ####
# Bring in contact rates
contact.est.slow.decay <- readRDS("Example simulation/Inputs/contact_rates_annual_lone_oaks.RDS")
contact.est.fast.decay <- readRDS("Example simulation/Inputs/contact_rates_annual_CP.RDS")

slow_kernel_df <-  tibble(distance_km = seq(0.001, 5, by = 0.001), 
                          contact = predict(contact.est.slow.decay, 
                                            newdata = tibble(distance = seq(0.001, 5, by = 0.001)), 
                                            re.form = NA, type = "response")) %>% 
  filter(contact >= 0.0001)

slow_contact_decay_sigma <- sigma(contact.est.slow.decay)

fast_kernel_df <- tibble(distance_km = seq(0.001, 3, by = 0.001), 
                         contact = predict(contact.est.fast.decay, 
                                           newdata = tibble(distance = seq(0.001, 3, by = 0.001)), 
                                           re.form = NA, type = "response")) %>% 
  filter(contact >= 0.0001)

fast_contact_decay_sigma <- sigma(contact.est.fast.decay)

# Set probability of infection given contact
sigma_value <- 0.018

# Set number of iterations for each study area x density x contact x season 
# combination. Right now, only one individual is infected at the beginning of
# each simulation.
starting.values <- set_simulation_parameters(n.sims = 1000, 
                                             init.prev = 0, 
                                             init.I = c(rep(1,1000)), 
                                             init.group = c(rep(0,1000)))


# Simulations ####
## Fall ####
### Elm Creek, MN (ECRP) ####
#### Low density ####
ECRP_fall_sim_input$pop_size <- (5)*sum(expanse(ECRP.grid, unit = "km")) #low density

RNGkind("L'Ecuyer-CMRG") #Sets random seed generator
set.seed(123)

##### Slow contact decay ####
initial_only_batch_result_lowD_slow_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                 sim.input = ECRP_fall_sim_input, 
                                                                 sigma_infection_prob = sigma_value, 
                                                                 grid = ECRP.grid, 
                                                                 centroids_df = ECRP_centroids_fall_df, 
                                                                 mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                 cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                 contact_distance_df = slow_kernel_df, 
                                                                 contact_distance_sigma = slow_contact_decay_sigma, 
                                                                 median_mvmt = ECRP_median_mvmt, 
                                                                 stop_I_0 = TRUE, 
                                                                 prop_residence_contact = 0, 
                                                                 rate_residence_contact = 0, 
                                                                 num_cores = 2)


##### Fast contact decay ####
initial_only_batch_result_lowD_fast_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                 sim.input = ECRP_fall_sim_input, 
                                                                 sigma_infection_prob = sigma_value, 
                                                                 grid = ECRP.grid, 
                                                                 centroids_df = ECRP_centroids_fall_df, 
                                                                 mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                 cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                 contact_distance_df = fast_kernel_df, 
                                                                 contact_distance_sigma = fast_contact_decay_sigma, 
                                                                 median_mvmt = ECRP_median_mvmt, 
                                                                 stop_I_0 = TRUE, 
                                                                 prop_residence_contact = 0, 
                                                                 rate_residence_contact = 0, 
                                                                 num_cores = 2)

#### Medium density ####
ECRP_fall_sim_input$pop_size <- (20)*expanse(ECRP_study_area, unit = "km") #low density

RNGkind("L'Ecuyer-CMRG") #Sets random seed generator
set.seed(123)

##### Slow contact decay ####
initial_only_batch_result_medD_slow_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                 sim.input = ECRP_fall_sim_input, 
                                                                 sigma_infection_prob = sigma_value, 
                                                                 grid = ECRP.grid, 
                                                                 centroids_df = ECRP_centroids_fall_df, 
                                                                 mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                 cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                 contact_distance_df = slow_kernel_df, 
                                                                 contact_distance_sigma = slow_contact_decay_sigma, 
                                                                 median_mvmt = ECRP_median_mvmt, 
                                                                 stop_I_0 = TRUE, 
                                                                 prop_residence_contact = 0, 
                                                                 rate_residence_contact = 0, 
                                                                 num_cores = 2)

##### Fast contact decay ####
initial_only_batch_result_medD_fast_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                 sim.input = ECRP_fall_sim_input, 
                                                                 sigma_infection_prob = sigma_value, 
                                                                 grid = ECRP.grid, 
                                                                 centroids_df = ECRP_centroids_fall_df, 
                                                                 mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                 cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                 contact_distance_df = fast_kernel_df, 
                                                                 contact_distance_sigma = fast_contact_decay_sigma, 
                                                                 median_mvmt = ECRP_median_mvmt, 
                                                                 stop_I_0 = TRUE, 
                                                                 prop_residence_contact = 0, 
                                                                 rate_residence_contact = 0, 
                                                                 num_cores = 2)

#### High density ####
ECRP_fall_sim_input$pop_size <- (50)*expanse(ECRP_study_area, unit = "km") #low density

RNGkind("L'Ecuyer-CMRG") #Sets random seed generator
set.seed(123)

##### Slow contact decay ####
initial_only_batch_result_highD_slow_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                  sim.input = ECRP_fall_sim_input, 
                                                                  sigma_infection_prob = sigma_value, 
                                                                  grid = ECRP.grid, 
                                                                  centroids_df = ECRP_centroids_fall_df, 
                                                                  mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                  cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                  contact_distance_df = slow_kernel_df, 
                                                                  contact_distance_sigma = slow_contact_decay_sigma, 
                                                                  median_mvmt = ECRP_median_mvmt, 
                                                                  stop_I_0 = TRUE, 
                                                                  prop_residence_contact = 0, 
                                                                  rate_residence_contact = 0, 
                                                                  num_cores = 2)

##### Fast contact decay ####
initial_only_batch_result_highD_fast_decay <- run_all_simulations_local_density_est(starting.values = starting.values, 
                                                                  sim.input = ECRP_fall_sim_input, 
                                                                  sigma_infection_prob = sigma_value, 
                                                                  grid = ECRP.grid, 
                                                                  centroids_df = ECRP_centroids_fall_df, 
                                                                  mvmt_dist_shape = ECRP_mvmt_dist_shape, 
                                                                  cand_pixels = ECRP_fall_dist_btw_pixels, 
                                                                  contact_distance_df = fast_kernel_df, 
                                                                  contact_distance_sigma = fast_contact_decay_sigma, 
                                                                  median_mvmt = ECRP_median_mvmt, 
                                                                  stop_I_0 = TRUE, 
                                                                  prop_residence_contact = 0, 
                                                                  rate_residence_contact = 0, 
                                                                  num_cores = 2)

#End of script