#-###########################-#
# Load libraries  ----
#-###########################-#
library(lmerTest)
library(emmeans)
library(glue)
library(dplyr)


#-###########################-#
# Define paths and global variable ----
#-###########################-#

## define globals
project_dir <- "/Users/aburns/Codebook/Projects/tda_criticality"      # Change to personal project directory
results_dir <- glue("{project_dir}/results/OASIS_EIB/longitudinal")
if (!file.exists(results_dir)) {dir.create(results_dir, recursive = TRUE)}
data_path <- glue("{project_dir}/data/From Igor_OASIS_Dataset")
scripts_path <- glue("{project_dir}/scripts")
src_path <- glue("{project_dir}/src")
image_path <- glue("{project_dir}/images/OASIS_Longitudinal_Images")
man_image_path <- glue("{image_path}/manuscript_images")


## Load in source files (helper functions and and addons):
source(glue("{src_path}/src_Linear modelling workflow_support functions.R")) 
source(glue("{src_path}/src_OASIS_EIB_utils.R"))
df_globalEIR <- read.csv("/Users/aburns/Codebook/Projects/tda_criticality/papers_and_protocols/MANUSCRIPT_OASIS_Longitudinal/CommBio_Submission/Accepted Version/SourceData_Figure3_3wayint.csv")
df_DMNLimbic <- read.csv("/Users/aburns/Codebook/Projects/tda_criticality/papers_and_protocols/MANUSCRIPT_OASIS_Longitudinal/CommBio_Submission/Accepted Version/SourceData_Figure4_DMNLimbic.csv")
#### Note: OASIS data use agreeement forbids seconday data sharing. We are allowed to share MR_ID and subject_ID 
#### of which a reader can use to replicate from the original datasets with their own dat use agreement


#-###########################-#
# Global EIR analysis ----
#-###########################-#

# Define model for Whole Brain EIR
model_wholebrain <- lmer(
  EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID), 
  data = df_globalEIR
)
summary(model_wholebrain)

# Post-hoc analysis for Whole Brain EIR
emmeans(model_wholebrain, pairwise ~ years_since_first_scan | sex)

# Generate interaction plot
plot_3way <- plot_model(
  model_wholebrain,
  type = "pred",
  terms = c("years_since_first_scan", "sex", "E4_binary"),
  ci.lvl = 0.95
) + 
  labs(
    title = "3-Way Interaction: Whole Brain EIR",
    x = "Years Since First Scan",
    y = "Estimated EIR"
  ) +
  theme_minimal()

# Save the plot
ggsave("WholeBrain_EIR_3WayInteraction.pdf", plot = plot_3way, width = 6, height = 4)
print(plot_3way)





#-###########################-#
# DMN-Limbic analysis ----
#-###########################-#

# Identify functional network EIR columns
func_cols <- grep("^func__", names(df_DMNLimbic), value = TRUE)

# Fit models for DMN-Limbic in-between EIR
results_list <- list()
for (col in func_cols[grepl("in_between", func_cols)]) {
  model <- tryCatch({
    lmer(as.formula(paste0("`", col, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")),
         data = df_DMNLimbic)
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    res <- as.data.frame(summary(model)$coefficients)
    res$EIR_type <- "in_between"
    res$Network <- col
    results_list[[col]] <- res
  }
}

# Estimated marginal means for DMN-Limbic in-between
col_interest <- "func__in_between__23"
model_dmn <- lmer(
  as.formula(paste0("`", col_interest, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")),
  data = df_DMNLimbic
)
em_dmn <- emmeans(model_dmn, ~ sex)
summary(em_dmn)

# Save plot of EIR estimates
library(ggplot2)
ggplot(as.data.frame(em_dmn), aes(x = sex, y = emmean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(title = "Estimated EIR by Sex", x = "Sex", y = "Estimated EIR") +
  theme_minimal()

ggsave("DMN_Limbic_EIR_Sex.pdf", width = 5, height = 4)
