# ---------*----------*---------*---------*---------*---------*---------*---------#
# ------------------------------------------------------------------------------- #
#                               Orientation Exp 1b.                               #
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
pacman::p_load(tidyverse, emmeans, tidyr, dplyr, knitr, lme4, afex, lmerTest, ggdist)
options(knitr.kable.NA = '') # hide NA with knitr function
devtools::install_github("psyteachr/introdataviz")
# ggplot2, carData, writexl, psych, afex, readr 

## ---------------------------------------------------------- #
## 1. Load data ####
## ---------------------------------------------------------- #
getwd()
#setwd("./mainExp")
datapath <- '../totalData/'
DataFiles <- list.files(path = datapath, pattern = "*.csv", full.names = TRUE)
results_table <- read.csv("results_table_long.csv")

maxmin_DataL <- data.frame()
sub_useL <- data.frame()


# Loop through each participant's data file
for (k in 1:length(DataFiles)) {
  
  # Load participant data
  participant_data <- read.csv(DataFiles[k])
  participant_data <- subset(participant_data, select = c("participant","age", "sex", "trial_phase","block","trial",
                               "image_name","imageCondition","rotationAngle","startOrientation","response_orientation","response_time","confidenceLevel", "imagePosition", "responseType"))
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
  grating_data$use <- ifelse(grating_data$corrected_response < 30, 1, 0) ##1 = use; 0 = don't use
  
  # Filter responsetime less than 800ms
  participant_data$rtUSE <- ifelse(participant_data$response_time > 800, 1, 0) ##1 = use; 0 = don't use
  rt_use <- participant_data %>%
    group_by(block) %>%
    summarise(
      participant = first(participant),
      rtUSESum = sum(rtUSE, na.rm = TRUE)
    )%>%
    mutate(
      rtUSE = ifelse(rtUSESum < 32, 0, 1)    # Set blockUSE to 0 if useSum < 0, otherwise 1
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
  
  # add imageCondition column
  ##filtered_participant_data <- filtered_participant_data %>%
    ##mutate(imageCondition = ifelse(grepl("max\\.png$", image_name), "max",
                                 ##ifelse(grepl("min\\.png$", image_name), "min", NA)))
  
  
  # Join with results_table on imgName (matching image index)
  filtered_participant_data <- filtered_participant_data %>%
    left_join(
      results_table %>% select(imgIdx_name, Mean_photo_deg, Mean_vecLD_deg), ##photo = filter; vecLD = contour
      by = c("image_name" = "imgIdx_name")
    )
  

  # Compare responseOrientation 
  maxmin_comparison <- filtered_participant_data %>%
    select(participant, image_name, imageCondition, confidenceLevel, imagePosition, responseType, corrected_response, Mean_photo_deg, Mean_vecLD_deg) %>%
    mutate(photo_comparison = pmin(abs(corrected_response - Mean_photo_deg),
                                   180 - abs(corrected_response - Mean_photo_deg)),
           vecLD_comparison = pmin(abs(corrected_response - Mean_vecLD_deg),
                                   180 - abs(corrected_response - Mean_vecLD_deg)))
  #maxmin_comparison <- unique(maxmin_comparison)

  # maxmin_comparison <- filtered_participant_data %>%
  #   select(participant, image_name, imageCondition, corrected_response, Mean_photo_deg, Mean_vecLD_deg) %>%
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

#trialsPerBLock = 80 (excluding grating)
#blockn = 8
100-sum(table(maxmin_DataL$participant))/(nrow(sub_info)*80*8)*100

maxmin_DataL_long <- maxmin_DataL %>%
  pivot_longer(cols = c(photo_comparison, vecLD_comparison), 
               names_to = "measure", 
               values_to = "value")

participant_means <- maxmin_DataL_long %>%
  group_by(participant, imageCondition, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = 'drop')

image_response_count <- maxmin_DataL_long %>%
  group_by(image_name) %>%
  summarise(response_count = n())

# Calculate the condition-level means and standard error
condition_summary <- maxmin_DataL_long %>%
  group_by(imageCondition, measure, responseType) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            confidenceM = mean(confidenceLevel, na.rm = TRUE),
            confidenceSD = sd(confidenceLevel, na.rm = TRUE),
            .groups = "drop")
condition_summary

# Marginal means for Image Condition (collapsed across measure)
imgCondition_summary <- maxmin_DataL_long %>%
  group_by(imageCondition) %>%
  summarise(M = mean(value, na.rm = TRUE),
            SEM = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Marginal means for Method 
method_summary <- maxmin_DataL_long %>%
  group_by(measure) %>%
  summarise(M = mean(value, na.rm = TRUE),
            SEM = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Marginal means for Response Type
responseType_summary <- maxmin_DataL_long %>%
  group_by(responseType) %>%
  summarise(M = mean(value, na.rm = TRUE),
            SEM = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

## ---------------------------------------------------------- #
## 3. ANOVA ####
## ---------------------------------------------------------- #
# 2-way repeated-measures ANOVA
anova_model <- aov_ez(
  id = "participant",          # subject identifier
  dv = "value",                # dependent variable
  data = maxmin_DataL_long,
  within = c("imageCondition", "measure", "responseType"),  # repeated measures factors
  factorize = FALSE            # factors are already factor type
)

# View ANOVA table
anova_model

# Post-hoc pairwise comparisons (measure within each imgCondition)
emmeans(anova_model, pairwise ~ measure | imageCondition, adjust = "tukey")

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
maxmin_DataL_long$participant <- as.factor(maxmin_DataL_long$participant)
maxmin_DataL_long$imageCondition <- as.factor(maxmin_DataL_long$imageCondition)
maxmin_DataL_long$measure <- as.factor(maxmin_DataL_long$measure)
maxmin_DataL_long$responseType <- as.factor(maxmin_DataL_long$responseType)

mixed_model <- lmer(value ~ (imageCondition * measure * responseType) + (1 + imageCondition * measure * responseType | participant), data = maxmin_DataL_long, 
                    contrasts=list(imageCondition=contr.sum, measure=contr.sum, responseType=contr.sum))
mixed_model2 <- lmer(value ~ imageCondition * measure * responseType + (1 + measure | participant), data = maxmin_DataL_long, 
                     contrasts=list(imageCondition=contr.sum, measure=contr.sum, responseType=contr.sum))
mixed_model3 <- lmer(value ~ imageCondition * measure * responseType + (1 + responseType | participant), data = maxmin_DataL_long, 
                     contrasts=list(imageCondition=contr.sum, measure=contr.sum, responseType=contr.sum))

AIC(mixed_model2, mixed_model3)
BIC(mixed_model2, mixed_model3)

summary(mixed_model2)

emmeans(mixed_model, pairwise ~ responseType | imageCondition + measure, pbkrtest.limit = 26136)

# install.packages("effectsize")  # Only needed once
# library(effectsize)
# eta_squared(mixed_model, partial = TRUE)
# anova_model <- aov(value ~ imageCondition * measure, 
#                    data = maxmin_DataL_long)
# eta_squared(anova_model, partial = TRUE)

# Plot the data

ggplot(participant_means, aes(x = measure, y = value, fill = imageCondition)) +
  # Split-violin plot to show distribution by measure within each imageCondition
  introdataviz::geom_split_violin(alpha = 0.4, trim = FALSE, show.legend = FALSE) +
  geom_errorbar(data = condition_summary, aes(x = measure, y = mean_value, ymin = mean_value - se_value, ymax = mean_value + se_value, 
                                              group = imageCondition), 
                width = 0.2, position = position_dodge(.175)) +
  geom_point(data = condition_summary, aes(x = measure, y = mean_value, group = imageCondition), show.legend = FALSE, 
             position = position_dodge(0.175), size = 1.5) +
  geom_point(data = participant_means, aes(x = measure, y = value, color = imageCondition, group = imageCondition), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.6) +
  geom_line(data = condition_summary, aes(x = measure, y = mean_value, group = imageCondition), 
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

ggplot(participant_means, aes(x = measure, y = value, fill = imageCondition)) +
  # Split-violin plot to show distribution by measure within each imageCondition
  introdataviz::geom_split_violin(alpha = 0.4, trim = FALSE, show.legend = FALSE) +
  geom_errorbar(data = condition_summary, aes(x = measure, y = mean_value, ymin = mean_value - se_value, ymax = mean_value + se_value, 
                                              group = imageCondition), 
                width = 0.2, position = position_dodge(.175)) +
  geom_point(data = condition_summary, aes(x = measure, y = mean_value, group = imageCondition), show.legend = FALSE, 
             position = position_dodge(0.175), size = 1.5) +
  geom_point(data = participant_means, aes(x = measure, y = value, color = imageCondition, group = imageCondition), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.6) +
  geom_line(data = condition_summary, aes(x = measure, y = mean_value, group = imageCondition), 
            position = position_dodge(0.175), linewidth = 0.7) +
  labs(x = "Method", y = "Response Error (degrees)") +
  scale_x_discrete(name = "Method", labels = c("Filter", "Contour")) +
  #scale_y_continuous(breaks=seq(20,70,10)) +
  scale_colour_discrete(name = "Image Condition", labels = c("Max", "Min"))+
  theme_classic() +
  ylim(5, 70)+
  guides(fill=FALSE)+
  theme(legend.position="none")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  theme(legend.title = element_text(size = 12))+
  theme(legend.text = element_text(size = 14)) +
  facet_wrap(~ responseType)+
  theme(strip.text = element_blank())

## ---------------------------------------------------------- #
## 5. Confidence Level ####
## ---------------------------------------------------------- #
# confidence level
maxmin_DataL %>% summarise(confidenceM = mean(confidenceLevel, na.rm = TRUE),
                           confidenceSD = sd(confidenceLevel, na.rm = TRUE),
                           .groups = 'drop')
maxmin_DataL %>%  
  group_by(imageCondition, responseType) %>%
  summarise(confidenceM = mean(confidenceLevel, na.rm = TRUE),
            confidenceSD = sd(confidenceLevel, na.rm = TRUE))

maxmin_DataL$participant <- as.factor(maxmin_DataL$participant)
maxmin_DataL$imageCondition <- as.factor(maxmin_DataL$imageCondition)
maxmin_DataL$responseType <- as.factor(maxmin_DataL$responseType)

kruskal_result <- kruskal.test(confidenceLevel ~ imageCondition,  data = maxmin_DataL)
print(kruskal_result)
kruskal_result <- kruskal.test(confidenceLevel ~ responseType,  data = maxmin_DataL)
print(kruskal_result)

maxmin_DataL_confidence <- maxmin_DataL_long %>%
  mutate(confidenceCategory = case_when(
    confidenceLevel %in% c(1,2) ~ "Low",
    confidenceLevel %in% c(4,5) ~ "High",
    TRUE ~ NA_character_  # Assign NA to other values
  )) %>%
  filter(!is.na(confidenceCategory))

maxmin_DataL_confidence %>%
  group_by(confidenceCategory, measure) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            confidenceM = mean(confidenceLevel, na.rm = TRUE),
            .groups = "drop")

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
maxmin_DataL_confidence$responseType <- as.factor(maxmin_DataL_confidence$responseType)
maxmin_DataL_confidence$measure <- as.factor(maxmin_DataL_confidence$measure)

mixed_model_confidence1 <- lmer(value ~ confidenceCategory * measure * responseType + (1 + confidenceCategory * measure * responseType | participant), data = maxmin_DataL_confidence)
mixed_model_confidence2 <- lmer(value ~ confidenceCategory * measure * responseType + (1 + responseType | participant), data = maxmin_DataL_confidence)
mixed_model_confidence3 <- lmer(value ~ confidenceCategory * measure * responseType + (1 + measure | participant), data = maxmin_DataL_confidence)
mixed_model_confidence4 <- lmer(value ~ confidenceCategory * measure * responseType + (1 + confidenceCategory | participant), data = maxmin_DataL_confidence)

AIC(mixed_model_confidence2, mixed_model_confidence3, mixed_model_confidence4)
BIC(mixed_model_confidence2, mixed_model_confidence3, mixed_model_confidence4)

summary(mixed_model_confidence4)

emmeans(mixed_model_confidence4, pairwise ~ confidenceCategory | measure, pbkrtest.limit = 26136)



