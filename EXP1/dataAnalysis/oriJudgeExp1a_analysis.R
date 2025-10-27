# ---------*----------*---------*---------*---------*---------*---------*---------#
# ------------------------------------------------------------------------------- #
#                               Orientation Exp 1a                                #
#
#            By Seohee Han
# ------------------------------------------------------------------------------- #
# Date: 2024-09-14
# Environment: R Studio Cloud, Windows 10 / macOS Big Sur
# ------------------------------------------------------------------------------- #
# ---------*----------*---------*---------*---------*---------*---------*---------#

# ===== CONTENTS ================================================================ #
#
# 1. Load data
# 2. Summary
# 3. ANOVA
# 4. Linear Mixed model
# 5. Confidence Level
#
# =============================================================================== #

# get ready
rm(list=ls())
set.seed(4228) # for replication

# load packages 
pacman::p_load(tidyverse, emmeans, tidyr, dplyr, knitr, lme4, afex, lmerTest, ggdist, pwr)
options(knitr.kable.NA = '') # hide NA with knitr function
devtools::install_github("psyteachr/introdataviz")
# ggplot2, carData, writexl, psych, afex, readr 

# For a two-sided paired t-test (simplest approximation)
n <- 31 #48  # number of participants
power <- 0.8
alpha <- 0.05

# Compute detectable Cohen's d
pwr.t.test(n = n, d = NULL, sig.level = alpha, power = power, type = "paired")$d


## ---------------------------------------------------------- #
## 1. Load data ####
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
    select(participant, image_name, imgCondition, corrected_response, confidenceLevel, Mean_photo_deg, Mean_vecLD_deg) %>%
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

## ---------------------------------------------------------- #
## 2. Summary ####
## ---------------------------------------------------------- #
sub_info <- unique(subset(sub_useL, select = c("participant","subUSE")))
table(maxmin_DataL$participant)

100-sum(table(maxmin_DataL$participant))/(nrow(sub_info)*464)*100

maxmin_DataL_long <- maxmin_DataL %>%
  pivot_longer(cols = c(photo_comparison, vecLD_comparison), 
               names_to = "measure", 
               values_to = "value")
#write.csv(maxmin_DataL_long, "maxmin_DataL_long.csv")

participant_means <- maxmin_DataL_long %>%
  group_by(participant, imgCondition, measure) %>%
  summarise(value = mean(value, na.rm = TRUE),
            confidenceM = mean(confidenceLevel, na.rm = TRUE), 
            .groups = 'drop')

image_response_count <- maxmin_DataL_long %>%
  group_by(image_name) %>%
  summarise(response_count = n())

# Calculate the condition-level means and standard error
condition_summary <- maxmin_DataL_long %>%
  group_by(imgCondition, measure) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            confidenceM = mean(confidenceLevel, na.rm = TRUE),
            .groups = "drop")
condition_summary

# Marginal means for Image Condition (collapsed across measure)
imgCondition_summary <- maxmin_DataL_long %>%
  group_by(imgCondition) %>%
  summarise(M = mean(value, na.rm = TRUE),
            SEM = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Marginal means for Method (collapsed across imgCondition)
method_summary <- maxmin_DataL_long %>%
  group_by(measure) %>%
  summarise(M = mean(value, na.rm = TRUE),
            SEM = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

maxmin_DataL_long$participant <- as.factor(maxmin_DataL_long$participant)
maxmin_DataL_long$imgCondition <- as.factor(maxmin_DataL_long$imgCondition)
maxmin_DataL_long$measure <- as.factor(maxmin_DataL_long$measure)


## ---------------------------------------------------------- #
## 3. ANOVA ####
## ---------------------------------------------------------- #
# 2-way repeated-measures ANOVA
anova_model <- aov_ez(
  id = "participant",          # subject identifier
  dv = "value",                # dependent variable
  data = maxmin_DataL_long,
  within = c("imgCondition", "measure"),  # repeated measures factors
  factorize = FALSE            # factors are already factor type
)

# View ANOVA table
anova_model

# Post-hoc pairwise comparisons (measure within each imgCondition)
emmeans(anova_model, pairwise ~ measure | imgCondition, adjust = "tukey")


## ---------------------------------------------------------- #
## 4. Linear Mixed model ####
## ---------------------------------------------------------- #
# mixed_model
#lmer default: treatment coding(dummy coding) 
# - the first level of a factor is taken as the baseline
# - the lower-order effects (such as main effects) are estimated at the level of the baseline, therefore yielding simple effects rather than main effects

#change to effects coding 
#for lmer (intercept = grand mean, estimates of the main factors*2 = difference between two levels)
#for lmer (if encode levels as -0.5 and +0.5 -> estimates of the main factors = difference between two levels)
#for lmer (example: contrasts(TotalData$Triplet)<-contr.sum(2)/2 #encoding your levels as -0.5 and +0.5)

# mixed_model <- lmer(value ~ (imgCondition * measure) + (1 | participant), data = maxmin_DataL_long, 
#                     contrasts=list(imgCondition=contr.sum, measure=contr.sum))
mixed_model2 <- lmer(value ~ imgCondition * measure + (1 + measure | participant), data = maxmin_DataL_long, 
                     contrasts=list(imgCondition=contr.sum, measure=contr.sum))
# mixed_model3 <- lmer(value ~ imgCondition * measure + (1 + imgCondition | participant), data = maxmin_DataL_long, 
#                      contrasts=list(imgCondition=contr.sum, measure=contr.sum))
# mixed_model4 <- lmer(value ~ imgCondition * measure + (1 + imgCondition *measure | participant), data = maxmin_DataL_long, 
#                      contrasts=list(imgCondition=contr.sum, measure=contr.sum))

# summary(mixed_model)
summary(mixed_model2)
# summary(mixed_model3)
# summary(mixed_model4)

emmeans(mixed_model2, pairwise ~ measure | imgCondition, pbkrtest.limit = 26136)

participant_means <- maxmin_DataL_long %>%
  group_by(participant, imgCondition, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

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
  ylim(5, 70)+
  theme_classic() +
  guides(fill=FALSE)+
  theme(legend.position = "right")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))+
  theme(legend.text = element_text(size = 12))

## ---------------------------------------------------------- #
## 5. Confidence Level ####
## ---------------------------------------------------------- #

# confidence level
maxmin_DataL %>% summarise(confidenceM = mean(confidenceLevel, na.rm = TRUE),
                           confidenceSD = sd(confidenceLevel, na.rm = TRUE),
                           .groups = 'drop')
maxmin_DataL %>%  
  group_by(imgCondition) %>%
  summarise(confidenceM = mean(confidenceLevel, na.rm = TRUE),
            confidenceSD = sd(confidenceLevel, na.rm = TRUE))
library(rstatix)   # For Kruskal-Wallis and Dunn's test

kruskal_result <- kruskal.test(confidenceLevel ~ imgCondition, data = maxmin_DataL)
print(kruskal_result)

maxmin_DataL_confidence <- maxmin_DataL_long %>%
  mutate(confidenceCategory = case_when(
    confidenceLevel %in% c(1,2) ~ "Low",
    confidenceLevel %in% c(4,5) ~ "High",
    TRUE ~ NA_character_  # Assign NA to other values
  )) %>%
  filter(!is.na(confidenceCategory))

maxmin_DataL_confidence %>%
  group_by(confidenceCategory) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            confidenceM = mean(confidenceLevel, na.rm = TRUE),
            .groups = "drop")

maxmin_DataL_confidence %>%
  group_by(confidenceCategory) %>%
  summarise(n = n())


maxmin_DataL_confidence$participant <- as.factor(maxmin_DataL_confidence$participant)
maxmin_DataL_confidence$confidenceCategory <- as.factor(maxmin_DataL_confidence$confidenceCategory)
maxmin_DataL_confidence$measure <- as.factor(maxmin_DataL_confidence$measure)

mixed_model_confidence1 <- lmer(value ~ confidenceCategory * measure + (measure | participant), data = maxmin_DataL_confidence)
mixed_model_confidence2 <- lmer(value ~ confidenceCategory * measure + (confidenceCategory | participant), data = maxmin_DataL_confidence)

AIC(mixed_model_confidence1, mixed_model_confidence2)
BIC(mixed_model_confidence1, mixed_model_confidence2)


summary(mixed_model_confidence2)
