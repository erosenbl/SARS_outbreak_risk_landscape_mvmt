#### Random Forest ####
# Load packages
library(ranger)
library(tibble)

# Read in data file
df_RF <- readRDS("Random forest analysis/Data/df_RF_final.RDS")

# Standardize landscape metrics by study area
df_RF <- df_RF %>%
group_by(study_area) %>%
  mutate(
    mean_origin_cell =
      ecdf(mean_origin_cell)(mean_origin_cell),
    
    sd_deer_origin_cell =
      ecdf(sd_deer_origin_cell)(sd_deer_origin_cell),
    
    mean_origin_kernel =
      ecdf(mean_origin_kernel)(mean_origin_kernel),
    
    sd_origin_kernel =
      ecdf(sd_origin_kernel)(sd_origin_kernel),
    
    sd_deer_origin_kernel =
      ecdf(sd_deer_origin_kernel)(sd_deer_origin_kernel),
    
    mean_origin_landscape =
      (mean_origin_landscape - mean(mean_origin_landscape)) / sd(mean_origin_landscape),
    
    sd_origin_landscape =
      (sd_origin_landscape - mean(sd_origin_landscape)) / sd(sd_origin_landscape),
    
    sd_deer_origin_landscape = (sd_deer_origin_landscape - mean(sd_deer_origin_landscape)) / sd(sd_deer_origin_landscape)) %>% 
  ungroup()


library(ranger)
library(DALEX)
  
covariates <- tibble(variable = c("rel_daily_SL_Male", "cv_sl_Female", "cv_sl_Male",
                                  "mean_origin_cell", "sd_deer_origin_cell", "frag_cell",
                                  "mean_origin_kernel", "sd_origin_kernel", "sd_deer_origin_kernel", "frag_kernel",
                                  "mean_origin_landscape", "sd_origin_landscape", "sd_deer_origin_landscape", "frag_landscape",
                                  "density_category", "contact_category"),
                     var_name = c("Male:female\ndisplacement ratio",
                                  "Female displacement\nheterogeneity",
                                  "Male displacement\nheterogeneity",
                                  "Mean centroid\nattractiveness",
                                  "Among-deer centroid\nheterogeneity",
                                  "Habitat fragmentation",
                                  "Mean centroid\nattractiveness",
                                  "Spatial heterogeneity of\ncentroid attractiveness",
                                  "Among-deer centroid\nheterogeneity",
                                  "Habitat fragmentation",
                                  "Mean centroid\nattractiveness",
                                  "Spatial heterogeneity of\ncentroid attractiveness",
                                  "Among-deer centroid\nheterogeneity",
                                  "Habitat fragmentation",
                                  "Deer density",
                                  "Contact kernel\ndecay rate"
                     ),
                     grouping = c("Animal movement", 
                                  "Animal movement", 
                                  "Animal movement",
                                  "Infection origin",
                                  "Infection origin",
                                  "Infection origin",
                                  "Index case's contact kernel",
                                  "Index case's contact kernel",
                                  "Index case's contact kernel",
                                  "Index case's contact kernel",
                                  "Study area",
                                  "Study area",
                                  "Study area", 
                                  "Study area",
                                  "Deer density and contact",
                                  "Deer density and contact"
                     ),
                     fill_col = c("#1B4F72","#1B4F72","#1B4F72",
                                  "#4A4A4A","#4A4A4A","#4A4A4A",
                                  "#287211","#287211","#287211","#287211",
                                  "#E3B778","#E3B778","#E3B778","#E3B778",
                                  "black", "black"
                     ))

covariates <- covariates %>% 
  mutate(grouping = factor(grouping, levels = c("Infection origin", "Index case's contact kernel", "Study area","Animal movement", "Deer density and contact")))

group_palette <- covariates %>%
  distinct(grouping, fill_col) %>%
  deframe()

##### Response Variable 1: Outbreak probability ####
outbreak_df_RF <- df_RF %>% 
  mutate(outbreak = factor(case_when(incidence_category != "Outbreak" ~ FALSE,
                                     TRUE ~ TRUE), levels = c(FALSE, TRUE)),
         rel_daily_SL_Male = daily_SL_Male/daily_SL_Female,
         frag_cell = 100-agg_index_cell,
         frag_kernel = 100-agg_index_kernel,
         frag_landscape = 100-agg_index_landscape) %>% 
  mutate(include_tag = case_when(density_category == "low" & contact_category == "slow" ~ 0,
                                 density_category == "low" & contact_category == "fast" ~ 0,
                                 TRUE ~ 1)) %>% 
  filter(include_tag == 1) %>% 
  dplyr::select(-include_tag)

variable <- "outbreak"

full_outbreak_df <- outbreak_df_RF %>% 
  mutate(density_contact = case_when(density_category == "low" & contact_category == "slow" ~ "Low density,\nlarge contact kernel",
                                     density_category == "low" & contact_category == "fast" ~ "Low density,\nsmall contact kernel",
                                     density_category == "med" & contact_category == "slow" ~ "Medium density,\nlarge contact kernel",
                                     density_category == "med" & contact_category == "fast" ~ "Medium density,\nsmall contact kernel",
                                     density_category == "high" & contact_category == "slow" ~ "High density,\nlarge contact kernel",
                                     density_category == "high" & contact_category == "fast" ~ "High density,\nlarge contact kernel")) %>% 
  mutate(density_category = case_when(density_category == "low" ~ "Low",
                                      density_category == "med" ~ "Medium",
                                      density_category == "high" ~ "High"),
         contact_category = case_when(contact_category == "slow" ~ "Large contact kernel",
                                      contact_category == "fast" ~ "Small contact kernel"),
         study_area = factor(study_area)) %>% 
  dplyr::select(variable,
                density_category,
                contact_category,
                rel_daily_SL_Male, 
                cv_sl_Female, 
                cv_sl_Male,
                mean_origin_cell, 
                sd_deer_origin_cell, 
                frag_cell,
                mean_origin_kernel, 
                sd_origin_kernel, 
                sd_deer_origin_kernel, 
                frag_kernel,
                mean_origin_landscape, 
                sd_origin_landscape, 
                sd_deer_origin_landscape, 
                frag_landscape,
                study_area)

n_features <- length(setdiff(names(full_outbreak_df), variable))

rf_full_outbreak_final <- ranger(
  outbreak ~ .,
  data = full_outbreak_df,
  mtry = 2,
  importance = "permutation",
  scale.permutation.importance = F,
  respect.unordered.factors = "order",
  min.node.size = 10,
  num.trees = 2000,
  sample.fraction = 0.5,
  replace = TRUE,
  seed = 123,
  probability = T
)

y_bin <- as.numeric(full_outbreak_df$outbreak == TRUE)

full_outbreak_explainer_ranger <- explain(
  model = rf_full_outbreak_final,
  data = full_outbreak_df[-1] %>% mutate(density_contact = case_when(density_category == "Low" & contact_category == "Large contact kernel" ~ "Low density,\nlarge contact kernel",
                                                                   density_category == "Low" & contact_category == "Small contact kernel" ~ "Low density,\nsmall contact kernel",
                                                                   density_category == "Medium" & contact_category == "Large contact kernel" ~ "Medium density,\nlarge contact kernel",
                                                                   density_category == "Medium" & contact_category == "Small contact kernel" ~ "Medium density,\nsmall contact kernel",
                                                                   density_category == "High" & contact_category == "Large contact kernel" ~ "High density,\nlarge contact kernel",
                                                                   density_category == "High" & contact_category == "Small contact kernel" ~ "High density,\nsmall contact kernel")),
  y = y_bin,
  label = "Ranger Model",
  verbose = FALSE
)

library(ggplot2)

# Importance plot
full_outbreak_importance_plot <- model_parts(full_outbreak_explainer_ranger, type = "difference", conditional = T) %>%
  as.data.frame(.) %>%
  group_by(variable) %>%
  summarise(
    mean_dropout_loss = mean(dropout_loss),
    .groups = "drop"
  ) %>%
  filter(!(variable %in% c("_baseline_", "_full_model_"))) %>% 
  merge(.,covariates) %>% 
  ggplot(., aes(reorder(var_name, mean_dropout_loss), mean_dropout_loss, fill = grouping)) +
  # geom_bar(stat="identity") +
  geom_col(position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Conditional permutation importance") +
  scale_fill_manual(
    values = group_palette,
    name = "Category"
  ) +
  coord_flip()+
  theme_classic()+
  labs(title = "Probability of outbreak")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text = element_text(size = 7),
        title = element_text(size = 8),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Effect plots for continous variables, using composite grouping of density and contact kernel

full_outbreak_ale_ranger <- model_profile(full_outbreak_explainer_ranger, type = "accumulated",
                                        groups = c("density_contact"), center = F)


# Save conditional effect plots for consolidation
full_outbreak_effect <- tibble(full_outbreak_ale_ranger$agr_profiles)

##### Response Variable 2: Maximum incidence rate ####
epi_df_RF <- df_RF %>% 
  filter(incidence_category == "Outbreak") %>% 
  mutate(max_incidence_rate_per10k = round(max_local_epi_rate*10000),
         rel_daily_SL_Male = daily_SL_Male/daily_SL_Female,
         frag_cell = 100-agg_index_cell,
         frag_kernel = 100-agg_index_kernel,
         frag_landscape = 100-agg_index_landscape) %>% 
  mutate(include_tag = case_when(density_category == "low" & contact_category == "slow" ~ 0,
                                 density_category == "low" & contact_category == "fast" ~ 0,
                                 TRUE ~ 1)) %>%
  filter(include_tag == 1) %>%
  dplyr::select(-include_tag)

variable <- "max_incidence_rate_per10k"

full_epi_df <- epi_df_RF %>% 
  mutate(density_contact = case_when(density_category == "low" & contact_category == "slow" ~ "Low density,\nlarge contact kernel",
                                     density_category == "low" & contact_category == "fast" ~ "Low density,\nsmall contact kernel",
                                     density_category == "med" & contact_category == "slow" ~ "Medium density,\nlarge contact kernel",
                                     density_category == "med" & contact_category == "fast" ~ "Medium density,\nsmall contact kernel",
                                     density_category == "high" & contact_category == "slow" ~ "High density,\nlarge contact kernel",
                                     density_category == "high" & contact_category == "fast" ~ "High density,\nlarge contact kernel")) %>% 
  mutate(density_category = case_when(density_category == "low" ~ "Low",
                                      density_category == "med" ~ "Medium",
                                      density_category == "high" ~ "High"),
         contact_category = case_when(contact_category == "slow" ~ "Large contact kernel",
                                      contact_category == "fast" ~ "Small contact kernel"),
         study_area = factor(study_area)) %>% 
  dplyr::select(variable,
                density_category,
                contact_category,
                rel_daily_SL_Male, 
                cv_sl_Female, 
                cv_sl_Male,
                mean_origin_cell, 
                sd_deer_origin_cell, 
                frag_cell,
                mean_origin_kernel, 
                sd_origin_kernel, 
                sd_deer_origin_kernel, 
                frag_kernel,
                mean_origin_landscape, 
                sd_origin_landscape, 
                sd_deer_origin_landscape, 
                frag_landscape,
                study_area)

n_features <- length(setdiff(names(full_epi_df), variable))

rf_full_epi_final <- ranger(
  log(max_incidence_rate_per10k) ~ .,
  data = full_epi_df,
  mtry = 3,
  importance = "permutation",
  scale.permutation.importance = F,
  respect.unordered.factors = "order",
  min.node.size = 15,
  num.trees = 2000,
  sample.fraction = 0.5,
  replace = TRUE,
  seed = 123
)

full_max_incidence_prop_rate_explainer_ranger <- explain(rf_full_epi_final, data = full_epi_df[-1] %>% mutate(density_contact = case_when(density_category == "Low" & contact_category == "Large contact kernel" ~ "Low density,\nlarge contact kernel",
                                                                                                                                     density_category == "Low" & contact_category == "Small contact kernel" ~ "Low density,\nsmall contact kernel",
                                                                                                                                     density_category == "Medium" & contact_category == "Large contact kernel" ~ "Medium density,\nlarge contact kernel",
                                                                                                                                     density_category == "Medium" & contact_category == "Small contact kernel" ~ "Medium density,\nsmall contact kernel",
                                                                                                                                     density_category == "High" & contact_category == "Large contact kernel" ~ "High density,\nlarge contact kernel",
                                                                                                                                     density_category == "High" & contact_category == "Small contact kernel" ~ "High density,\nsmall contact kernel")), y = log(full_epi_df$max_incidence_rate_per10k), label = "Ranger Model", verbose = FALSE)

# Importance plot
full_max_incidence_prop_rate_importance_plot <- model_parts(full_max_incidence_prop_rate_explainer_ranger, type = "difference", conditional = T) %>%
  as.data.frame(.) %>%
  group_by(variable) %>%
  reframe(
    mean_dropout_loss = mean(dropout_loss),
    .groups = "drop"
  ) %>%
  filter(!(variable %in% c("_baseline_", "_full_model_"))) %>%
  merge(.,covariates) %>%
  ggplot(., aes(reorder(var_name, mean_dropout_loss), mean_dropout_loss, fill = grouping)) +
  # geom_bar(stat="identity") +
  geom_col(position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Conditional permutation importance") +
  scale_fill_manual(
    values = group_palette,
    name = "Category"
  ) +
  coord_flip()+
  labs(title = "Maximum incidence rate")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text = element_text(size = 7),
        title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

full_epi_ale_ranger <- model_profile(full_max_incidence_prop_rate_explainer_ranger, type = "accumulated", center = F, groups = c("density_contact"))

# Save conditional effect plots for consolidation
full_epi_effect <- tibble(full_epi_ale_ranger$agr_profiles)

##### Response Variable 3: Maximum spatial spread ####
spatial_df_RF <- df_RF %>% 
  filter(incidence_category == "Outbreak") %>% 
  mutate(rel_daily_SL_Male = daily_SL_Male/daily_SL_Female,
         frag_cell = 100-agg_index_cell,
         frag_kernel = 100-agg_index_kernel,
         frag_landscape = 100-agg_index_landscape) %>% 
  mutate(include_tag = case_when(density_category == "low" & contact_category == "slow" ~ 0,
                                 density_category == "low" & contact_category == "fast" ~ 0,
                                 #density_category == "med" & contact_category == "fast" ~ 0,
                                 TRUE ~ 1)) %>% 
  filter(include_tag == 1) %>% 
  dplyr::select(-include_tag)

variable <- "max_spatial_spread"

full_spatial_df <- spatial_df_RF %>% 
  mutate(density_contact = case_when(density_category == "low" & contact_category == "slow" ~ "Low density,\nlarge contact kernel",
                                     density_category == "low" & contact_category == "fast" ~ "Low density,\nsmall contact kernel",
                                     density_category == "med" & contact_category == "slow" ~ "Medium density,\nlarge contact kernel",
                                     density_category == "med" & contact_category == "fast" ~ "Medium density,\nsmall contact kernel",
                                     density_category == "high" & contact_category == "slow" ~ "High density,\nlarge contact kernel",
                                     density_category == "high" & contact_category == "fast" ~ "High density,\nlarge contact kernel")) %>% 
  mutate(density_category = case_when(density_category == "low" ~ "Low",
                                      density_category == "med" ~ "Medium",
                                      density_category == "high" ~ "High"),
         contact_category = case_when(contact_category == "slow" ~ "Large contact kernel",
                                      contact_category == "fast" ~ "Small contact kernel"),
         study_area = factor(study_area)) %>% 
  dplyr::select(variable,
                density_category,
                contact_category,
                rel_daily_SL_Male, 
                cv_sl_Female, 
                cv_sl_Male,
                mean_origin_cell, 
                sd_deer_origin_cell, 
                frag_cell,
                mean_origin_kernel, 
                sd_origin_kernel, 
                sd_deer_origin_kernel, 
                frag_kernel,
                mean_origin_landscape, 
                sd_origin_landscape, 
                sd_deer_origin_landscape, 
                frag_landscape,
                study_area)


rf_full_spatial_final <- ranger(
  log(max_spatial_spread+0.001) ~ .,
  data = full_spatial_df,
  mtry = 2,
  importance = "permutation",
  scale.permutation.importance = F,
  respect.unordered.factors = "order",
  min.node.size = 10,
  num.trees = 2000,
  sample.fraction = 0.5,
  replace = TRUE,
  seed = 123
)

full_max_spatial_rate_explainer_ranger <- explain(rf_full_spatial_final, data = full_spatial_df[-1] %>% mutate(density_contact = case_when(density_category == "Low" & contact_category == "Large contact kernel" ~ "Low density,\nlarge contact kernel",
                                                                                                                                           density_category == "Low" & contact_category == "Small contact kernel" ~ "Low density,\nsmall contact kernel",
                                                                                                                                           density_category == "Medium" & contact_category == "Large contact kernel" ~ "Medium density,\nlarge contact kernel",
                                                                                                                                           density_category == "Medium" & contact_category == "Small contact kernel" ~ "Medium density,\nsmall contact kernel",
                                                                                                                                           density_category == "High" & contact_category == "Large contact kernel" ~ "High density,\nlarge contact kernel",
                                                                                                                                           density_category == "High" & contact_category == "Small contact kernel" ~ "High density,\nsmall contact kernel")), 
                                                  y = log(full_spatial_df$max_spatial_spread+0.001), label = "Ranger Model", verbose = FALSE)

# Importance plot
full_max_spatial_rate_importance_plot <- model_parts(full_max_spatial_rate_explainer_ranger, type = "difference", conditional = T) %>%
  as.data.frame(.) %>%
  group_by(variable) %>%
  reframe(
    mean_dropout_loss = mean(dropout_loss),
    .groups = "drop"
  ) %>%
  filter(!(variable %in% c("_baseline_", "_full_model_"))) %>%
  merge(.,covariates) %>%
  ggplot(., aes(reorder(var_name, mean_dropout_loss), mean_dropout_loss, fill = grouping)) +
  # geom_bar(stat="identity") +
  geom_col(position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Conditional permutation importance") +
  scale_fill_manual(
    values = group_palette,
    name = "Category"
  ) +
  coord_flip()+
  labs(title = "Maximum spatial spread rate")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text = element_text(size = 7),
        title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

full_spatial_ale_ranger <- model_profile(full_max_spatial_rate_explainer_ranger, type = "accumulated",
                                     groups = c("density_contact"), center = F)

# Save conditional effect plots for consolidation
full_spatial_effect <- tibble(full_spatial_ale_ranger$agr_profiles)
