# ---------*----------*---------*---------*---------*---------*---------*---------#
# ------------------------------------------------------------------------------- #
#                               R2                          #
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
pacman::p_load(tidyverse, emmeans, tidyr, dplyr, knitr, ggpubr, rstatix, lmerTest, ggdist)
options(knitr.kable.NA = '') # hide NA with knitr function

## ---------------------------------------------------------- #
## 1. load data ####
## ---------------------------------------------------------- #
getwd()
#setwd("./mainExp")
# datapath <- '../data'
# DataFiles <- list.files(path = datapath, pattern = "*.csv", full.names = TRUE)
R2_old <- read.csv("allroiR2old.csv", header = FALSE)
R2_control <- read.csv("allroiR2control.csv", header = FALSE)
R2_new <- read.csv("allroiR2ori.csv", header = FALSE)
R2_control_sfmean <- read.csv("allroiR2_sfmean_control.csv", header = FALSE)
R2_old_sfmean <- read.csv("allroiR2_sfmean_old.csv", header = FALSE)
# colnames(R2_old) <- c("old")  # Add the appropriate column names
# colnames(R2_control) <- c("control")  # Add the appropriate column names
colnames(R2_new) <- c("new")  # Add the appropriate column names
colnames(R2_control_sfmean) <- c("control")
colnames(R2_old_sfmean) <- c("old")

allR2 <-cbind(R2_old_sfmean, R2_control_sfmean, R2_new)
controlR2 <-cbind(R2_control, R2_control_sfmean)

allR2_long <- allR2 %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "R2")
controlR2_long <- controlR2 %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "R2")

kruskal_result <- kruskal.test(R2 ~ Method, data = allR2_long)
kruskal_result
kruskal_result_control <- kruskal.test(R2 ~ Method, data = controlR2_long)
kruskal_result_control

pairwise_results <- list()
# Loop through all pairwise combinations of conditions
conditions <- unique(allR2_long$Method)
for (i in 1:(length(conditions)-1)) {
  for (j in (i+1):length(conditions)) {
    # Extract the values for each condition
    group1 <- allR2_long$R2[allR2_long$Method == conditions[i]]
    group2 <- allR2_long$R2[allR2_long$Method == conditions[j]]
    
    # Perform Wilcoxon test for each pair and store the W-statistic and p-value
    test_result <- wilcox.test(group1, group2)
    pairwise_results[[paste(conditions[i], conditions[j], sep = "_vs_")]] <- test_result
  }
}

# View the W-statistics for each pairwise comparison
for (comparison in names(pairwise_results)) {
  cat("Comparison:", comparison, "\n")
  cat("W-statistic:", pairwise_results[[comparison]]$statistic, "\n")
  cat("p-value:", pairwise_results[[comparison]]$p.value, "\n\n")
}


#plot 
condition_summary <- allR2_long %>%
group_by(Method) %>%
  summarise(mean_value = mean(R2, na.rm = TRUE),
            sd_value = sd(R2, na.rm = TRUE),
            se_value = sd(R2, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
condition_summary$Method <- factor(condition_summary$Method, 
                                   levels = c("old", "control", "new"))

ggplot(allR2_long, aes(x = Method, y = R2, fill = Method))+
  #geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1, outliers = FALSE)+
  theme_classic()



# Bar plot with custom labels and colors
ggplot(condition_summary, aes(x = Method, y = mean_value, fill = Method)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2) +
  scale_x_discrete(labels = c("old" = "Photo\nSteerable Pyramid", 
                              "control" = "Line drawing\nSteerable Pyramid", 
                              "new" = "Contour")) +
  scale_fill_manual(values = c("old" = "#0070C0", 
                               "control" = "#4EA72E", 
                               "new" = "#E54291")) +
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position="none",
    axis.text = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.title = element_text(size = 14, family = "Helvetica", color = "black")
  )+
  labs(x = NULL, y = NULL)

filtered_summary <- condition_summary %>%
  filter(Method %in% c("old", "new"))

# Bar plot for photo-steerable pyramid and contour
ggplot(filtered_summary, aes(x = Method, y = mean_value, fill = Method)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2) +
  labs(x = "Method",
       y = "Mean R2") +
  scale_x_discrete(labels = c("old" = "photo\nsteerable pyramid", 
                              "new" = "contour")) +
  scale_fill_manual(values = c("old" = "#0070C0", 
                               "new" = "#E54291")) +
  theme_classic() +
  theme(legend.position="none")+
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))

#old 0070C0
#control 4EA72E
#contour E54291

subMeanR2_rhcontrol <- read.csv("subMeanR2_rhcontrol.csv", header = FALSE)
subMeanR2_lhcontrol <- read.csv("subMeanR2_lhcontrol.csv", header = FALSE)
subMeanR2_rhold <- read.csv("subMeanR2_rhold.csv", header = FALSE)
subMeanR2_lhold <- read.csv("subMeanR2_lhold.csv", header = FALSE)
subMeanR2_rhori <- read.csv("subMeanR2_rhori.csv", header = FALSE)
subMeanR2_lhori <- read.csv("subMeanR2_lhori.csv", header = FALSE)

subMeanR2_control <- rbind(subMeanR2_rhcontrol, subMeanR2_lhcontrol)
subMeanR2_old <- rbind(subMeanR2_rhold, subMeanR2_lhold)
subMeanR2_ori <- rbind(subMeanR2_rhori, subMeanR2_lhori)

colnames(subMeanR2_old) <- c("old")  # Add the appropriate column names
colnames(subMeanR2_control) <- c("control")  # Add the appropriate column names
colnames(subMeanR2_ori) <- c("new")  # Add the appropriate column names

allsubMeanR2 <-cbind(subMeanR2_old, subMeanR2_control, subMeanR2_ori)

allsumbMeanR2_long <- allsubMeanR2 %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "R2")

# Plots
ggplot(allsumbMeanR2_long, aes(x = Method, y = R2, fill = Method))+
#  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1, outliers = FALSE)+
  theme_classic()




