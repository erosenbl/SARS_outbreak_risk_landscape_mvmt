# Script to designate spatial contact kernels as slow or fast decay, and 
# identify one slow decay and fast decay kernels with the same area under the
# curve. 

# This script is for demonstration purposes only - contact rate files are not 
# included in this repo.

# Function code
library(glmmTMB)
library(DescTools)
library(dplyr)

# Annual contact rates
# MN
ECRP_annual <- readRDS("~contact_rates_annual_EC.RDS")
CP_annual <- readRDS("~contact_rates_annual_CP.RDS")
shakopee_annual <- readRDS("~contact_rates_annual_shakopee.RDS")

# IL
TON_annual <- readRDS("~contact_rates_annual_TON.RDS")
Shelbyville_annual <- readRDS("~contact_rates_annual_Shelby.RDS")

# TN
ames_annual <- readRDS("~contact_rates_annual_ames.RDS")
lone_oaks_annual <- readRDS("~contact_rates_annual_lone_oaks.RDS")

# NY
SI_annual <- readRDS("~contact_rates_annual_SI.RDS")

# PA
TL_annual <- readRDS("~contact_rates_annual_TL.RDS")

annual_models <- list(ECRP_annual = ECRP_annual, 
                      CP_annual = CP_annual, 
                      shakopee_annual = shakopee_annual, 
                      TON_annual = TON_annual, 
                      Shelbyville_annual = Shelbyville_annual, 
                      ames_annual = ames_annual, 
                      lone_oaks_annual = lone_oaks_annual, 
                      SI_annual = SI_annual, 
                      TL_annual = TL_annual)

# standard distance range
distance_seq <- seq(0, 2, by = 0.01)

# Calculate pairwise contact rate for each model
predict_contact <- function(model, dist_seq) {
  newdata <- data.frame(distance = dist_seq)
  preds <- predict(model, newdata = newdata, type = "response", re.form = NA)
  data.frame(distance = dist_seq, pred = preds)
}

pred_list <- map(annual_models, predict_contact, dist_seq = distance_seq)

calc_half_life <- function(pred_df) {
  max_pred <- max(pred_df$pred, na.rm = TRUE)
  half_val <- max_pred / 2
  # find smallest distance where pred <= half_val
  idx <- which(pred_df$pred <= half_val)[1]
  if (is.na(idx)) return(NA)
  pred_df$distance[idx]
}

half_lives <- map_dbl(pred_list, calc_half_life)

#Classify decay rate as fast or slow, relative to median half life across models
decay_category <- ifelse(half_lives >= median(half_lives, na.rm = TRUE), "slow",
                         "fast")

#Calculate area under curve
auc_values <- map_dbl(pred_list, ~ AUC(x = .$distance, y = .$pred))

# Combine metrics into table
model_stats <- tibble(
  model = names(annual_models),
  half_life = half_lives,
  auc = auc_values,
  decay = decay_category
) %>% 
  sort(.,"auc")

# Explore pairings
slow_models <- model_stats %>% filter(decay == "slow")
fast_models <- model_stats %>% filter(decay == "fast")

pairings <- map_dfr(slow_models$model, function(sm) {
  slow_auc <- slow_models$auc[slow_models$model == sm]
  fast_match <- fast_models$model[which.min(abs(fast_models$auc - slow_auc))]
  tibble(slow_model = sm,
         fast_model = fast_match,
         auc_diff = abs(fast_models$auc[fast_models$model == fast_match] - 
                          slow_auc))
})

pairings


library(ggplot2)

pred_df <- bind_rows(
  map2(pred_list, names(pred_list), ~ mutate(.x, model = .y))
) %>%
  left_join(model_stats, by = "model")


ggplot(pred_df, aes(x = distance, y = pred, color = decay, group = model)) +
  geom_line(alpha = 0.6) +
  labs(x = "Distance", y = "Predicted Mean", color = "Decay Type") +
  theme_minimal()

pred_df_top_pair <- pred_df %>% 
  filter(model %in% c("lone_oaks_annual", "CP_annual"))

ggplot(pred_df, aes(x = distance, y = pred, color = decay, group = model)) +
  geom_line(alpha = 0.3) +
  geom_line(data = pred_df_top_pair, aes(x = distance, y = pred, color = decay, 
                                         group = model),
            linewidth = 1.3) +
  labs(x = "Distance", y = "Predicted Mean", color = "Decay Type") +
  scale_colour_viridis_d(end = 0.7)+
  theme_classic()


# Define the sequence of distances

# Predict function that handles zero-inflated Gamma models with random effects
predict_glmmTMB_ci <- function(model, dist_seq, include_re = FALSE) {
  newdata <- data.frame(distance = dist_seq)
  
  # Predict on the link scale to get correct SEs
  preds <- predict(model, newdata = newdata, type = "response", 
                   re.form = if (include_re) NULL else NA)
  
  # Back-transform and compute confidence intervals
  tibble(
    distance = dist_seq,
    est = preds,
    shape = 1/(sigma(model)^2),
    scale = est/shape,
    pred_lower = qgamma(p = 0.025, shape = shape, scale = scale),
    pred_upper = qgamma(p = 0.975, shape = shape, scale = scale)
  )
}

# Apply to all models in your list
pred_list <- map(annual_models, predict_glmmTMB_ci, dist_seq = distance_seq, 
                 include_re = FALSE)

# Combine into one tidy data frame
pred_df <- bind_rows(
  map2(pred_list, names(pred_list), ~ mutate(.x, model = .y))
)

# (Optional) join decay info or metadata
pred_df <- left_join(pred_df, model_stats, by = "model")

highlight_models <- c("lone_oaks_annual", "CP_annual")

pred_df_top_pair <- pred_df %>% 
  filter(model %in% highlight_models)

contact_rate_all_plot <- ggplot(pred_df, 
                                aes(x = distance, y = est, color = decay, 
                                    group = model)) +
  geom_line(alpha = 0.6) +
  labs(x = "Distance (km)", 
       y = "Pair-wise contacts per day",
       color = "Decay Type", 
       fill = "Decay Type") +
  scale_colour_viridis_d(end = 0.7) +
  scale_fill_viridis_d(end = 0.7) +
  scale_x_continuous(
    limits = c(-1.25, max(pred_df$distance)),
    breaks = c(0, 0.5, 1.0, 1.5, 2.0)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1.0)
  ) +
  coord_cartesian(clip = "off") +   # allow drawing outside panel
  theme_classic() +
  theme(
    plot.margin = margin(5.5, 30, 5.5, 60)  # extra left margin
  )

label_df <- pred_df %>%
  group_by(model) %>%
  filter(distance == min(distance)) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(study_area = c("Carver Park (MN)", 
                        "Elm Creek Reserve\nPark (MN)",
                        "Staten Island (NY)",
                        "Shelbyville (IL)",
                        "Tresaure Lake (PA)",
                        "Touch of Nature (IL)",
                        "Ames (TN)",
                        "Lone Oaks (TN)",
                        "Shakopee (MN)")) %>% 
  arrange(desc(est)) %>%
  mutate(y_lab = c(1,
                   0.95,
                   0.55,
                   0.47,
                   0.42,
                   0.37,
                   0.34,
                   0.3,
                   0.23))

contact_rate_all_plot <- contact_rate_all_plot +
  geom_text(
    data = label_df,
    aes(
      x = -0.02,     # position left of curves
      y = y_lab,
      label = study_area,
      color = decay
    ),
    hjust = 1,
    size = 3,
    show.legend = FALSE
  )

top_pair_label <- label_df %>% filter(study_area %in% c("Carver Park (MN)", 
                                                        "Lone Oaks (TN)"))

contact_rate_highlight_plot <- ggplot(pred_df, aes(x = distance, y = est, 
                                                   color = decay, group = model)
                                      ) +
  geom_line(alpha = 0.6) +
  geom_ribbon(data = pred_df_top_pair,
              aes(ymin = pred_lower, ymax = pred_upper, fill = decay),
              alpha = 0.25, color = NA) +
  geom_line(data = pred_df_top_pair,
            aes(x = distance, y = est, color = decay, group = model),
            linewidth = 1.3) +
  labs(x = "Distance (km)", y = "Pair-wise contacts per day", 
       color = "Decay Type", fill = "Decay Type") +
  geom_text(
    data = top_pair_label,
    aes(
      x = -0.02,     # position left of curves
      y = est,
      label = study_area,
      color = decay
    ),
    hjust = 1,
    size = 3,
    show.legend = FALSE
  )+
  scale_colour_viridis_d(end = 0.7)+
  scale_fill_viridis_d(end = 0.7)+
  scale_x_continuous(
      limits = c(-0.75, max(pred_df$distance)),
      breaks = c(0, 0.5, 1.0, 1.5, 2.0)
    ) +
  scale_y_continuous(
      breaks = c(0, 0.5, 1),
      limits = c(0, 1.25)
    ) +
  theme_classic()
  
#Going to go with lone oaks and carver for slow and fast, respectively
