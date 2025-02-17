#=============================================================================
#
# Evaluating effects of data quality and variable weighting on habitat suitability modelling 
# Script 3) Habitat suitability index (HSI) modelling
#
# Authors: Stephanie Arsenault, Robyn Linner, and Yong Chen
# Date: 02.12.2025
#
## LOAD REQUIRED PACKAGES =====================================================
#
l_packages <- c("stats", "gbm", "tidyr", "dplyr", "ggplot2", "mgcv", "tidyverse", "sf", "colorRamps")
for (p in l_packages){
  if(! p %in% installed.packages()){
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = T)
}
#

## Step 0- DATA PREP  =========================================================
# 
# Load data
uncalib_cpue <- read.csv() # add path - uncalibrated CPUE
calib_cpue <- read.csv() # add path - calibrated CPUE

calib.weight <- read.csv() # add path - annual HSI values generated using the calibrated weighted model 
calib.noweight <- read.csv() # add path - annual HSI values generated using the calibrated unweighted model 

uncalib.weight <- read.csv() # add path - annual HSI values generated using the uncalibrated weighted model 
uncalib.noweight <- read.csv() # add path - annual HSI values generated using the uncalibrated unweighted model 


# Process function
process_data <- function(data, value_name) {
  data %>%
    pivot_longer(c(3:32), names_to = "Year_Month", values_to = value_name) %>%
    separate(Year_Month, into = c("Month", "Year"), sep = "_") %>%
    select(-Month)
}

# Process uncalibrated data 
uncalib.weight <- process_data(uncalib.weight, "Uncalib.weight")
uncalib.noweight <- process_data(uncalib.noweight, "Uncalib.noweight")
uncalib.models <- left_join(uncalib.weight, uncalib.noweight, by = c("lat", "lon", "Year"))
uncalib.models.cpue <- left_join(uncalib.models, uncalib_cpue, by = c("lat", "lon", "Year")) %>% na.omit()

# Process calibrated data
calib.weight <- process_data(calib.weight, "Calib.weight")
calib.noweight <- process_data(calib.noweight, "Calib.noweight")
calib.models <- left_join(calib.noweight, calib.weight, by = c("lat", "lon", "Year"))
calib.models.cpue <- left_join(calib.models, calib_cpue, by = c("lat", "lon", "Year")) %>% na.omit()
#
## Step 1- CROSS VALIDATION  ==================================================
# 
# Cross-validation function
cross_validate <- function(data, formula, num_iterations = 100) {
  metrics <- data.frame(Iteration = 1:num_iterations, RMSE = numeric(num_iterations), AIC = numeric(num_iterations),
                        RRMSE = numeric(num_iterations), R2 = numeric(num_iterations), slope = numeric(num_iterations),
                        intercept = numeric(num_iterations))
  
  for (i in 1:num_iterations) {
    train_indices <- sample(seq_len(nrow(data)), size = 0.8 * nrow(data))
    test_indices <- setdiff(seq_len(nrow(data)), train_indices)
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    model <- lm(formula, data = train_data)
    predicted_values <- predict(model, newdata = test_data)
    metrics[i, c("intercept", "slope")] <- coef(model)
    metrics[i, "R2"] <- summary(model)$r.squared
    metrics[i, "RMSE"] <- sqrt(mean((test_data[[all.vars(formula)[1]]] - predicted_values)^2))
    metrics[i, "AIC"] <- AIC(model)
    metrics[i, "RRMSE"] <- metrics[i, "RMSE"] / mean(test_data[[all.vars(formula)[1]]])
  }
  
  metrics
}

# Uncalibrated models
metrics_uncal_noweight <- cross_validate(uncalib.models.cpue, Uncalib.noweight ~ cpue)
metrics_uncal_weight <- cross_validate(uncalib.models.cpue, Uncalib.weight ~ cpue)

# Calibrated models
metrics_cal_noweight <- cross_validate(calib.models.cpue, Calib.noweight ~ cal.data)
metrics_cal_weight <- cross_validate(calib.models.cpue, Calib.weight ~ cal.data)

# Save results
write.csv(metrics_uncal_noweight, "metrics_uncal_noweight.csv")
write.csv(metrics_uncal_weight, "metrics_uncal_weight.csv")
write.csv(metrics_cal_noweight, "metrics_cal_noweight.csv")
write.csv(metrics_cal_weight, "metrics_cal_weight.csv")


# Mean metrics and 95% confidence intervals

# Load metrics files
metrics_files <- list(
  "metrics_uncal_noweight.csv",
  "metrics_uncal_weight.csv",
  "metrics_cal_noweight.csv",
  "metrics_cal_weight.csv"
)

# Read data
dataframes <- lapply(metrics_files, read.csv)

# Function to calculate mean and 95% confidence interval
calc_stats <- function(df) {
  means <- colMeans(df, na.rm = TRUE)
  ci_lower <- apply(df, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  ci_upper <- apply(df, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  data.frame(means, ci_lower, ci_upper)
}

# Apply function to each dataframe and combine results
stats <- do.call(rbind, lapply(dataframes, calc_stats))

# Save results
write.csv()
#
## Step 1- ANNUAL PERCETAGES OF POOR, FAIR, GOOD HABITAT  =====================
#
# Load data
calib.weight <- read.csv() # add path - annual HSI values generated using the calibrated weighted model 
calib.noweight <- read.csv() # add path - annual HSI values generated using the calibrated unweighted model 

uncalib.weight <- read.csv() # add path - annual HSI values generated using the uncalibrated weighted model 
uncalib.noweight <- read.csv() # add path - annual HSI values generated using the uncalibrated unweighted model 

# Function to process data
process_data <- function(data) {
  data %>%
    pivot_longer(cols = 3:32, names_to = "Year_Month", values_to = "HSI") %>%
    separate(Year_Month, into = c("Month", "Year"), sep = "_") %>%
    select(-Month)
}

# Process data
calib.weight <- process_data(calib.weight)
calib.noweight <- process_data(calib.noweight)
uncalib.weight <- process_data(uncalib.weight)
uncalib.noweight <- process_data(uncalib.noweight)

# Function to calculate proportions
calculate_proportions <- function(data, label) {
  data %>%
    group_by(Year) %>%
    summarise(
      proportion_less_than_0_3 = sum(HSI < 0.3, na.rm = TRUE) / sum(!is.na(HSI)),
      proportion_between_0_3_and_0_7 = sum(HSI >= 0.3 & HSI < 0.7, na.rm = TRUE) / sum(!is.na(HSI)),
      proportion_greater_than_or_equal_to_0_7 = sum(HSI >= 0.7, na.rm = TRUE) / sum(!is.na(HSI)),
      .groups = 'drop'
    ) %>%
    drop_na() %>%
    mutate(Data = label)
}

# Calculate proportions
proportion_calib_weight <- calculate_proportions(calib.weight, "Calibrated-weighted")
proportion_calib_noweight <- calculate_proportions(calib.noweight, "Calibrated-unweighted")
proportion_uncalib_weight <- calculate_proportions(uncalib.weight, "Uncalibrated-weighted")
proportion_uncalib_noweight <- calculate_proportions(uncalib.noweight, "Uncalibrated-unweighted")

# Combine the data frames and reshape to longer format
proportions <- bind_rows(proportion_calib_weight, proportion_calib_noweight, proportion_uncalib_weight, proportion_uncalib_noweight) %>%
  pivot_longer(
    cols = starts_with("proportion"),
    names_to = "proportion_type",
    values_to = "proportion"
  )

proportions$Year <- as.numeric(proportions$Year)


# PLOTS
# Convert 'proportion_type' to a factor and specify the levels
proportions$proportion_type <- factor(proportions$proportion_type, 
                                      levels = c("proportion_greater_than_or_equal_to_0_7",
                                                "proportion_between_0_3_and_0_7",
                                                "proportion_less_than_0_3"))

proportions <- proportions %>%
  mutate(proportion_type = recode(proportion_type,
                                  "proportion_greater_than_or_equal_to_0_7" = "Good", 
                                  "proportion_between_0_3_and_0_7" = "Fair", 
                                  "proportion_less_than_0_3"= "Poor"))

# Plot
proportions$proportion_type <- as.factor(proportions$proportion_type)
proportions$Data <- factor(proportions$Data, levels = c("Calibrated-weighted", "Calibrated-unweighted", 
                                                        "Uncalibrated-weighted", "Uncalibrated-unweighted"))

proportions_summary <- proportions %>%
  group_by(Data, proportion_type) %>%
  summarise(mean_proportion = mean(proportion), se_proportion = sd(proportion) / sqrt(n()))

# Perform linear regression for each group and calculate p-values
lm_results <- proportions %>%
  group_by(Data, proportion_type) %>%
  do({
    model <- lm(proportion ~ Year, data = .)
    model_summary <- summary(model)
    data.frame(p_value = model_summary$coefficients[2, 4], r_squared = model_summary$r.squared)
  })

# Add a linetype column based on the significance of the regression
lm_results <- lm_results %>%
  mutate(linetype = ifelse(p_value < 0.05, "dashed", "solid"))

# Merge the linetype information back to the original data
proportions <- proportions %>%
  left_join(lm_results %>% select(Data, proportion_type, linetype), by = c("Data", "proportion_type"))

my_labeller <- function(variable, value) {
  return(paste0(LETTERS[1:length(value)], ") ", value))
}

# Plot all data points and add linear trendlines with confidence intervals for significant relationships
ggplot(proportions, aes(x = Year, y = proportion, color = proportion_type, group = proportion_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(linetype = linetype, color = proportion_type, fill = proportion_type), size = 0.5, show.legend = FALSE) +
  facet_wrap(~ Data, nrow = 2, labeller = my_labeller) +
  scale_color_manual(values = c("Good" = "#377F37", "Fair" = "#7CD1FF", "Poor" = "#FFD700")) +
  scale_fill_manual(values = c("Good" = "#377F37", "Fair" = "#7CD1FF", "Poor" = "#FFD700")) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
  labs(x = "Year", y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.5, "cm"), panel.spacing.y = unit(0.2, "cm"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(family = "Arial", color = "black", size = 12),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(color = "black", family = "Arial", size = 12),
        strip.text = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        text = element_text(family = "Arial"),
        strip.background = element_blank(),
        legend.margin = margin(t = -8),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  scale_y_continuous(expand = c(0.03, 0.03)) +
  scale_x_continuous(expand = c(0.03, 0.03))
#

## Step 2- MEDIAN HSI COMPARISON PLOT =========================================
#Load data
calib.weight <- read.csv() # add path - annual HSI values generated using the calibrated weighted model 
calib.noweight <- read.csv() # add path - annual HSI values generated using the calibrated unweighted model 

uncalib.weight <- read.csv() # add path - annual HSI values generated using the uncalibrated weighted model 
uncalib.noweight <- read.csv() # add path - annual HSI values generated using the uncalibrated unweighted model 

# Function to calculate median HSI and remove columns
process_median <- function(data) {
  data %>%
    rowwise() %>%
    mutate(median_HSI = median(c_across(3:32), na.rm = TRUE)) %>%
    ungroup() %>%
    select(-c(3:32))
}

# Process data
calib.weight <- process_median(calib.weight) %>% mutate(Data = "Calibrated- weighted")
calib.noweight <- process_median(calib.noweight) %>% mutate(Data = "Calibrated- unweighted")
uncalib.weight <- process_median(uncalib.weight) %>% mutate(Data = "Uncalibrated- weighted")
uncalib.noweight <- process_median(uncalib.noweight) %>% mutate(Data = "Uncalibrated- unweighted")

# Combine data
HSI.med <- bind_rows(calib.noweight, calib.weight, uncalib.noweight, uncalib.weight)

# Set factor levels
HSI.med$Data <- factor(HSI.med$Data, levels = c("Calibrated- weighted", "Calibrated- unweighted", 
                                                "Uncalibrated- weighted", "Uncalibrated- unweighted"))


midpoint.x <- mean(c(-74.05, -73.85))
midpoint.y<-mean(c(40.75,41.3))

labeller_func <- function(variable,value){
  return(paste0(str_replace_all(value, " ", "\n")))
}

combined_labeller <- function(variable, value) {
  labels <- paste0(LETTERS[1:length(value)], ") ", str_replace_all(value, " ", "\n"))
  return(labels)
}

(hsi_facet = ggplot(data = HSI.med,aes(x=lon,y=lat))+
    geom_tile(aes(fill = median_HSI))+
    scale_fill_gradientn(colours = matlab.like(100), limits = c(0, 1), "") +
    coord_quickmap(xlim = c(-74.05,-73.85), ylim = c(40.7,41.35)) +
    labs(x = "Longitude (°W)", y = "Latitude (°N)")+
    scale_x_continuous(expand = c(0,0), breaks = c(-74, -73.9)) +
    scale_y_continuous(expand = c(0,0), breaks = c(40.75, 41.03,41.3))+ 
    facet_wrap(~Data, nrow=1, ncol=4, labeller = combined_labeller)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 12, family = "Arial"),
          strip.background = element_blank(),
          axis.text.x = element_text(family = "Arial", color="black", size = 10),
          axis.text.y = element_text(color="black", family="Arial", size=10),
          axis.title.x = element_text(size = 12, family = "Arial"),
          axis.title.y = element_text(size = 12, family = "Arial"),
          panel.spacing.x = unit(0.5, "lines"),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

