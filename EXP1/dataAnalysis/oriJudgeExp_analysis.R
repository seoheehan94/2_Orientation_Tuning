# ---------*----------*---------*---------*---------*---------*---------*---------#
# ------------------------------------------------------------------------------- #
#                               Orientation Exp -pilot                            #
#
#            By Seohee Han
# ------------------------------------------------------------------------------- #
# Date: 2024-09-14
# Environment: R Studio Cloud, Windows 10 / macOS Big Sur
# ------------------------------------------------------------------------------- #
# ---------*----------*---------*---------*---------*---------*---------*---------#

# -- REQUIRED PACKAGES -- #
# tidyverse, emmeans, knitr, psych, dplyr, afex, BayesFactor,
# car, ggplot2, lme4, lmerTest, simr, sjPlot, tibble, readr, pwr2, sjPlot, effects
# -- Download packages before start -- #


# ===== CONTENTS ================================================================ #
#
# 1. load data
# 2. Combine subjects' data
# 3. Filter blocks under cutoff accuracy
# 4. Summary
# 5. Accuracy Analysis
#   5.1. one sample t-test
#   5.2. Wilcoxon rank test
#   5.3. Generalized Linear Mixed Modeling 
# 6. Memorability analysis
#   6.1. Consistency between participants
#   6.2. LD vs Photo memorability

# =============================================================================== #

# get ready
rm(list=ls())
set.seed(4228) # for replication

# load packages 
pacman::p_load(tidyverse, emmeans, tidyr, dplyr, knitr, lme4, afex, lmerTest, ggdist)
options(knitr.kable.NA = '') # hide NA with knitr function
devtools::install_github("psyteachr/introdataviz")
# ggplot2, carData, writexl, psych, afex, readr 

## ---------------------------------------------------------- #
## 1. load data ####
## ---------------------------------------------------------- #
getwd()
#setwd("./mainExp")
datapath <- '../data'
DataFiles <- list.files(path = datapath, pattern = "*.csv", full.names = TRUE)
results_table <- read.csv("../rawStimuli/results_table_long.csv")

maxmin_DataL <- data.frame()
sub_useL <- data.frame()


# Loop through each participant's data file
for (k in 1:length(DataFiles)) {
  
  # Load participant data
  participant_data <- read.csv(DataFiles[k])
  participant_data <- subset(participant_data, select = c("participant","trial_phase","block","trial",
                               "image_name","rotationAngle","startOrientation","response_orientation","bar_response_time","confidenceLevel"))
  participant_data <- participant_data %>% filter(trial_phase == 'main')
  participant_data$image_name <- gsub('stimuli/', '', participant_data$image_name)
 
  # Filter catch trials
  grating_data <- participant_data %>%
    filter(grepl("grating", image_name))
  grating_data <- unique(grating_data)
  grating_data <- grating_data %>%
    mutate(rotationAngle2 = rotationAngle %% 180)
  grating_data <- grating_data %>%
    mutate(corrected_response =  abs(response_orientation - rotationAngle2)) %>%
    mutate(corrected_response = ifelse(corrected_response > 90, abs(corrected_response - 180), corrected_response))
  grating_data$use <- ifelse(grating_data$corrected_response < 30, 1, 0)
  
  # Filter responsetime less than 800ms
  participant_data$rtUSE <- ifelse(participant_data$bar_response_time > 800, 1, 0)
  rt_use <- participant_data %>%
    group_by(block) %>%
    summarise(
      participant = first(participant),
      rtUSESum = sum(rtUSE, na.rm = TRUE)
    )%>%
    mutate(
      rtUSE = ifelse(rtUSESum < 29, 0, 1)    # Set blockUSE to 0 if useSum < 0, otherwise 1
    ) %>%
    ungroup()
  
  grating_use <- grating_data %>%
    group_by(block) %>%
    summarise(
      participant = first(participant),      # Take the first value of 'participant' in each group
      useSum = sum(use, na.rm = TRUE)        # Summarize 'use' within each 'block'
    ) %>%
    mutate(
      gratingUSE = ifelse(useSum < 1, 0, 1)    # Set blockUSE to 0 if useSum < 1, otherwise 1
    ) %>%
    ungroup() 
  
  sub_use <- rt_use %>%
    select(block, participant, rtUSE) %>%
    left_join(
      grating_use %>% select(block, participant, gratingUSE), 
      by = c("block", "participant")
    ) %>%
    mutate(
      blockUSE = ifelse(rtUSE == 1 & gratingUSE == 1, 1, 0)  # Correct condition for blockUSE
    ) %>%
    group_by(participant) %>%  # Group by participant to calculate sum of blockUSE for each participant
    mutate(
      subUSE = ifelse(sum(blockUSE) < 5, 0, 1)  # Calculate subUSE based on the sum of blockUSE
    ) %>%
    ungroup()
  
  sub_useL <- rbind(sub_useL, sub_use)
  participant_data <- participant_data %>%
    left_join(sub_use %>% select(participant, block, blockUSE, subUSE), by = c("participant", "block"))
  filtered_participant_data <- participant_data %>%
    filter(blockUSE == 1 & subUSE == 1 & rtUSE ==1)
  filtered_participant_data <- filtered_participant_data %>%
    filter(!grepl("grating", image_name))
  
  # Adjust responseOrientation based on rotationAngle
  filtered_participant_data <- filtered_participant_data %>%
    mutate(corrected_response = (response_orientation - rotationAngle) %% 360) %>%
    mutate(corrected_response = ifelse(corrected_response > 180, corrected_response - 180, corrected_response))
  
  # add imgCondition column
  filtered_participant_data <- filtered_participant_data %>%
    mutate(imgCondition = ifelse(grepl("max\\.png$", image_name), "max",
                                 ifelse(grepl("min\\.png$", image_name), "min", NA)))
  
  
  # Join with results_table on imgName (matching image index)
  filtered_participant_data <- filtered_participant_data %>%
    left_join(
      results_table %>% select(imgIdx_name, Mean_photo_deg, Mean_vecLD_deg), 
      by = c("image_name" = "imgIdx_name")
    )
  

  # Compare responseOrientation 
  maxmin_comparison <- filtered_participant_data %>%
    select(participant, image_name, imgCondition, corrected_response, Mean_photo_deg, Mean_vecLD_deg) %>%
    mutate(photo_comparison = pmin(abs(corrected_response - Mean_photo_deg),
                                   180 - abs(corrected_response - Mean_photo_deg)),
           vecLD_comparison = pmin(abs(corrected_response - Mean_vecLD_deg),
                                   180 - abs(corrected_response - Mean_vecLD_deg)))
  maxmin_comparison <- unique(maxmin_comparison)

  # maxmin_comparison <- filtered_participant_data %>%
  #   select(participant, image_name, imgCondition, corrected_response, Mean_photo_deg, Mean_vecLD_deg) %>%
  #   mutate(photo_comparison = ifelse(abs(corrected_response - Mean_photo_deg) <= abs(180 - abs(corrected_response - Mean_photo_deg)), 
  #                                    corrected_response - Mean_photo_deg, 
  #                                    (180 - abs(corrected_response - Mean_photo_deg)) * sign(corrected_response - Mean_photo_deg)),
  #          vecLD_comparison = ifelse(abs(corrected_response - Mean_vecLD_deg) <= abs(180 - abs(corrected_response - Mean_vecLD_deg)), 
  #                                    corrected_response - Mean_vecLD_deg, 
  #                                    (180 - abs(corrected_response - Mean_vecLD_deg)) * sign(corrected_response - Mean_vecLD_deg)))
  # 
  maxmin_DataL <- rbind(maxmin_DataL, maxmin_comparison)

}



# summary
sub_info <- unique(subset(sub_useL, select = c("participant","subUSE")))
table(maxmin_DataL$participant)

100-sum(table(maxmin_DataL$participant))/(nrow(sub_info)*464)*100

maxmin_DataL_long <- maxmin_DataL %>%
  pivot_longer(cols = c(photo_comparison, vecLD_comparison), 
               names_to = "measure", 
               values_to = "value")

participant_means <- maxmin_DataL_long %>%
  group_by(participant, imgCondition, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = 'drop')

image_response_count <- maxmin_DataL_long %>%
  group_by(image_name) %>%
  summarise(response_count = n())

# Calculate the condition-level means and standard error
condition_summary <- maxmin_DataL_long %>%
  group_by(imgCondition, measure) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
condition_summary



# mixed_model
#lmer default: treatment coding(dummy coding) 
# - the first level of a factor is taken as the baseline
# - the lower-order effects (such as main effects) are estimated at the level of the baseline, therefore yielding simple effects rather than main effects

#change to effects coding 
#for lmer (intercept = grand mean, estimates of the main factors*2 = difference between two levels)
#for lmer (if encode levels as -0.5 and +0.5 -> estimates of the main factors = difference between two levels)
#for lmer (example: contrasts(TotalData$Triplet)<-contr.sum(2)/2 #encoding your levels as -0.5 and +0.5)
maxmin_DataL_long$participant <- as.factor(maxmin_DataL_long$participant)
maxmin_DataL_long$imgCondition <- as.factor(maxmin_DataL_long$imgCondition)
maxmin_DataL_long$measure <- as.factor(maxmin_DataL_long$measure)

mixed_model <- lmer(value ~ (imgCondition * measure) + (1 | participant), data = maxmin_DataL_long, 
                    contrasts=list(imgCondition=contr.sum, measure=contr.sum))
mixed_model2 <- lmer(value ~ imgCondition * measure + (1 + measure | participant), data = maxmin_DataL_long, 
                     contrasts=list(imgCondition=contr.sum, measure=contr.sum))

summary(mixed_model)
summary(mixed_model2)

emmeans(mixed_model, pairwise ~ measure | imgCondition, pbkrtest.limit = 26136)


# Plot the data

ggplot(participant_means, aes(x = measure, y = value, fill = imgCondition)) +
  # Split-violin plot to show distribution by measure within each imgCondition
  introdataviz::geom_split_violin(alpha = 0.4, trim = FALSE, show.legend = FALSE) +
  geom_errorbar(data = condition_summary, aes(x = measure, y = mean_value, ymin = mean_value - se_value, ymax = mean_value + se_value, 
                                              group = imgCondition), 
                width = 0.2, position = position_dodge(.175)) +
  geom_point(data = condition_summary, aes(x = measure, y = mean_value, group = imgCondition), show.legend = FALSE, 
             position = position_dodge(0.175), size = 1.5) +
  geom_point(data = participant_means, aes(x = measure, y = value, color = imgCondition, group = imgCondition), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.6) +
  geom_line(data = condition_summary, aes(x = measure, y = mean_value, group = imgCondition), 
            position = position_dodge(0.175), linewidth = 0.7) +
  labs(x = "Method", y = "Response Error (degrees)") +
  scale_x_discrete(name = "Method", labels = c("Filter", "Contour")) +
  scale_colour_discrete(name = "Image Condition", labels = c("Max", "Min"))+
  ylim(13, 65)+
  theme_classic() +
  guides(fill=FALSE)+
  theme(legend.position = "right")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))+
  theme(legend.text = element_text(size = 12))


### why large range for min?
#High
# 12
# 24
# Sona25
# Sona32
# Sona69
# 
# Low
# 13
# Sona19
# Sona24
# Sona60
# Sona61
# Filter data for specific participants
high_data <- maxmin_DataL_long %>%
  filter(participant %in% c("12", "24", "sona25", "sona32", "sona69"))
low_data <- maxmin_DataL_long %>%
  filter(participant %in% c("13", "sona19", "sona24", "sona60", "sona61"))

high_summary <- high_data %>%
  group_by(imgCondition, measure) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
low_summary <- low_data %>%
  group_by(imgCondition, measure) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
high_summary
low_summary
# Filter data for specific participants, imageCondition "min", and measure "photo_comparison"
high_data_min <- maxmin_DataL_long %>%
  filter(participant %in% c("12", "24", "sona25", "sona32", "sona69"),
         imgCondition == "min",
         measure == "photo_comparison")%>%
  mutate(participant = as.factor(participant))
low_data_min <- maxmin_DataL_long %>%
  filter(participant %in% c("13", "sona19", "sona24", "sona60", "sona61"),
         imgCondition == "min",
         measure == "photo_comparison")%>%
  mutate(participant = as.factor(participant))

# Plot the distribution of 'value'
ggplot(high_data_min, aes(x = value, fill = participant)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of 'value' high_data_min",
       x = "Value",
       y = "Frequency",
       fill = "Participant") +
  theme_minimal() +
  ylim(0,45)+
  scale_fill_brewer(palette = "Set2") # Optional for better colors

ggplot(low_data_min, aes(x = value, fill = participant)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of 'value' low_data_min",
       x = "Value",
       y = "Frequency",
       fill = "Participant") +
  theme_minimal() +
  ylim(0,45)+
  scale_fill_brewer(palette = "Set2")
#
overlapping_names <- distinct(data.frame(intersect(high_data_min$image_name, low_data_min$image_name)))

# proportion of corrected_response
bin_width <- 10
high_data_over70 <- high_data_min %>%
  filter(value >= 70) %>%
  mutate(bin = cut(corrected_response, 
                   breaks = seq(0, 180, by = bin_width), 
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

high_data_less70 <- high_data_min %>%
  filter(value < 70) %>%
  mutate(bin = cut(corrected_response, 
                   breaks = seq(0, 180, by = bin_width), 
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

ggplot(high_data_over70, aes(x = bin, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "high_data_over70 - Corrected Response",
       x = "Corrected Response Bins",
       y = "Proportion") +
  theme_minimal() +
  ylim(0,0.35)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(high_data_less70, aes(x = bin, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "high_data_less70 - Corrected Response",
       x = "Corrected Response Bins",
       y = "Proportion") +
  theme_minimal() +
  ylim(0,0.35)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

low_data_over70 <- low_data_min %>%
  filter(value >= 70) %>%
  mutate(bin = cut(corrected_response, 
                   breaks = seq(0, 180, by = bin_width), 
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))
low_data_less70 <- low_data_min %>%
  filter(value < 70) %>%
  mutate(bin = cut(corrected_response, 
                   breaks = seq(0, 180, by = bin_width), 
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

ggplot(low_data_over70, aes(x = bin, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "low_data_over70 - Corrected Response",
       x = "Corrected Response Bins",
       y = "Proportion") +
  theme_minimal() +
  ylim(0,0.35)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(low_data_less70, aes(x = bin, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "low_data_over70 - Corrected Response",
       x = "Corrected Response Bins",
       y = "Proportion") +
  theme_minimal() +
  ylim(0,0.35)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


high_data_max <- maxmin_DataL_long %>%
  filter(participant %in% c("12", "24", "sona25", "sona32", "sona69"),
         imgCondition == "max",
         measure == "photo_comparison")%>%
  mutate(participant = as.factor(participant))
low_data_max <- maxmin_DataL_long %>%
  filter(participant %in% c("13", "sona19", "sona24", "sona60", "sona61"),
         imgCondition == "max",
         measure == "photo_comparison")%>%
  mutate(participant = as.factor(participant))

# Plot the distribution of 'value'
ggplot(high_data_max, aes(x = value, fill = participant)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of 'value' high_data_max",
       x = "Value",
       y = "Frequency",
       fill = "Participant") +
  theme_minimal() +
  ylim(0,45)+
  scale_fill_brewer(palette = "Set2") # Optional for better colors

ggplot(low_data_max, aes(x = value, fill = participant)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of 'value' low_data_max",
       x = "Value",
       y = "Frequency",
       fill = "Participant") +
  theme_minimal() +
  ylim(0,45)+
  scale_fill_brewer(palette = "Set2")

#img12408_min
#img12859_min
#img11369_min
#img14572_min

#makes sense
#img15938_min
#img13850_min
  
ggplot(maxmin_DataL_long, aes(x = Mean_photo_deg)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of Mean_photo_deg",
       x = "Value",
       y = "Frequency") +
  theme_minimal() 
ggplot(maxmin_DataL_long, aes(x = Mean_vecLD_deg)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "dodge") +
  labs(title = "Distribution of Mean_vecLD_deg",
       x = "Value",
       y = "Frequency") +
  theme_minimal() 

