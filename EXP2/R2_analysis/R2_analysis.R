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

R2_contour <- read.csv("allroiR2__contour.csv", header = FALSE)
R2_ldSP <- read.csv("allroiR2__ldSP.csv", header = FALSE)
R2_photoSP <- read.csv("allroiR2__photoSP.csv", header = FALSE)

R2_all <- bind_rows(
  mutate(R2_contour, Method = "contour", ROI = "all"),
  mutate(R2_ldSP, Method = "ldSP", ROI = "all"),
  mutate(R2_photoSP, Method = "photoSP", ROI = "all")
)
colnames(R2_all)[1] <- "R2"


versions <- paste0("V", 1:4)
types <- c("contour", "ldSP", "photoSP")
colnames_map <- c("contour", "ldSP", "photoSP")

for (v in versions) {
  for (i in seq_along(types)) {
    type <- types[i]
    name <- colnames_map[i]
    file <- paste0(v, "R2", type, ".csv")
    varname <- paste0("R2_", name, "_", v)
    
    data <- read.csv(file, header = FALSE)
    colnames(data) <- "R2" 
    
    data <- data %>%
      mutate(Method = name,
             ROI = paste0(v)  )
    
    assign(varname, data)
  }
}

version_dfs <- mget(ls(pattern = "^R2_.*_V[1-4]$"))
R2_ROI <- bind_rows(version_dfs)
R2_all_ROI <- bind_rows(R2_all, R2_ROI)


## ---------------------------------------------------------- #
## 2. Statistical test ####
## ---------------------------------------------------------- #

kruskal_result <- kruskal.test(R2 ~ Method, data = R2_all)
kruskal_result

pairwise_results <- list()

# Loop through all pairwise combinations of conditions
conditions <- unique(R2_all$Method)
for (i in 1:(length(conditions)-1)) {
  for (j in (i+1):length(conditions)) {
    # Extract the values for each condition
    group1 <- R2_all$R2[R2_all$Method == conditions[i]]
    group2 <- R2_all$R2[R2_all$Method == conditions[j]]
    
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

# Get unique ROIs
ROIs <- unique(R2_all_ROI$ROI)

# Initialize lists to store results
kruskal_results <- list()
pairwise_results <- list()

for (roi in ROIs) {
  cat("=== ROI:", roi, "===\n")
  
  # Subset data for this ROI
  data_roi <- subset(R2_all_ROI, ROI == roi)
  
  # Kruskal-Wallis test
  kruskal_results[[roi]] <- kruskal.test(R2 ~ Method, data = data_roi)
  print(kruskal_results[[roi]])
  
  # Pairwise Wilcoxon tests
  conditions <- unique(data_roi$Method)
  pairwise_results[[roi]] <- list()
  
  for (i in 1:(length(conditions)-1)) {
    for (j in (i+1):length(conditions)) {
      group1 <- data_roi$R2[data_roi$Method == conditions[i]]
      group2 <- data_roi$R2[data_roi$Method == conditions[j]]
      
      test_result <- wilcox.test(group1, group2)
      pairwise_results[[roi]][[paste(conditions[i], conditions[j], sep = "_vs_")]] <- test_result
    }
  }
  
  # Print pairwise results for this ROI
  for (comparison in names(pairwise_results[[roi]])) {
    cat("Comparison:", comparison, "\n")
    cat("W-statistic:", pairwise_results[[roi]][[comparison]]$statistic, "\n")
    cat("p-value:", pairwise_results[[roi]][[comparison]]$p.value, "\n\n")
  }
}

## ---------------------------------------------------------- #
## 3. Plots ####
## ---------------------------------------------------------- #
# Bar plot with custom labels and colors
condition_summary <- R2_all %>%
  group_by(Method) %>%
  summarise(mean_value = mean(R2, na.rm = TRUE),
            sd_value = sd(R2, na.rm = TRUE),
            se_value = sd(R2, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
condition_summary$Method <- factor(condition_summary$Method, 
                                   levels = c("photoSP", "ldSP", "contour"))

ggplot(condition_summary, aes(x = Method, y = mean_value, fill = Method)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2) +
  scale_x_discrete(labels = c("photoSP" = "Photo\nSteerable Pyramid", 
                              "ldSP" = "Line drawing\nSteerable Pyramid", 
                              "contour" = "Contour")) +
  scale_fill_manual(values = c("photoSP" = "#0070C0", 
                               "ldSP" = "#4EA72E", 
                               "contour" = "#E54291")) +
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position="none",
    axis.text = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.title = element_text(size = 14, family = "Helvetica", color = "black")
  )+
  labs(x = NULL, y = NULL)

# plot by ROI
condition_summary_ROI <- R2_all_ROI %>%
  group_by(ROI, Method) %>%
  summarise(
    mean_value = mean(R2, na.rm = TRUE),
    sd_value = sd(R2, na.rm = TRUE),
    se_value = sd(R2, na.rm = TRUE) / sqrt(n())
  )

condition_summary_ROI$Method <- factor(condition_summary_ROI$Method, 
                                   levels = c("photoSP", "ldSP", "contour"))
condition_summary_ROI$ROI <- factor(condition_summary_ROI$ROI, 
                                       levels = c("all", "V1", "V2", "V3", "V4"))
levels(condition_summary_ROI$ROI)[levels(condition_summary_ROI$ROI) == "V4"] <- "hV4"

ggplot(condition_summary_ROI, aes(x = Method, y = mean_value, fill = Method)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                width = 0.2) +
  facet_wrap(~ ROI, nrow = 1) +  # keep ROI titles
  scale_fill_manual(values = c("photoSP" = "#0070C0",
                               "ldSP" = "#4EA72E",
                               "contour" = "#E54291")) +
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),   # remove x-axis labels
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    axis.text.y = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.title = element_text(size = 14, family = "Helvetica", color = "black"),
    strip.background = element_blank(),   # remove gray facet box
    strip.text = element_text(size = 14, family = "Helvetica", color = "black")
  ) +
  labs(x = NULL, y = NULL)
