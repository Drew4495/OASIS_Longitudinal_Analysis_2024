#---------------------------------------------------#
#    LMMs on Functional Network EIR (DMN and Limbic aggregated)
#    Replicating workflow for EIR types: "in", "between", "in_between"
#---------------------------------------------------#



#-##############################################-#
# Define Globals and load in packages and src ----
#-##############################################-#
rm(list = ls()) 

## Load in libraries
library(glue)
library(glmmTMB)
library(dplyr)
library(nortest)
library(moments)
library(ggpubr)
library(lmerTest)
library(purrr)
library(sjPlot)
library(gridExtra)
library(rmarkdown)
library(knitr)
library(showtext)
library(openxlsx)
font_add("Arial", regular= "arial.ttf")
showtext_auto()


## define globals
project_dir <- "/Users/aburns/Codebook/Projects/tda_criticality"      # Change to personal project directory
results_dir <- glue("{project_dir}/results/OASIS_EIB")
if (!file.exists(results_dir)) {dir.create(results_dir, recursive = TRUE)}
data_path <- glue("{project_dir}/data/From Igor_OASIS_Dataset")
scripts_path <- glue("{project_dir}/scripts")
src_path <- glue("{project_dir}/src")
image_path <- glue("{project_dir}/images/OASIS_Longitudinal_Images")
man_image_path <- glue("{image_path}/manuscript_images")


## Load in source files (helper functions and and addons):
source(glue("{src_path}/src_Linear modelling workflow_support functions.R")) 
source(glue("{src_path}/src_OASIS_EIB_utils.R"))





#-###########################-#
# Just DMN: Import df and Preprocess  ----
#-###########################-#

## Import data frame 
oasis_df_path <- glue("{data_path}/preprocessed/oasis_df_all_apoenpscore_CDRbinary_dupes_allregions.csv")
df <- read.csv(oasis_df_path, header=T)[,-1]

## import other dataframes with functional info to merge with above df
filepath_func_df <- paste0(data_path, "/OASIS_Data_Func_Networks_DMNaggregated_EIR.xlsx")
func_df = read.xlsx(filepath_func_df)


### Merge func_df with demographic df
func_df$subject_idx = func_df$subject_idx + 1
colnames(df)[colnames(df) == "SUBJECT."] <- "subject_idx"
df <- merge(df, func_df, by="subject_idx", all.x = TRUE)

## Rename columns
df <- df %>%
  dplyr::rename(
    MR_ID = `MR.ID`,
    subject_ID = Subject,
    EIR_whole_brain = Whole.Brain,
    EIR_right_hippo = Right.Hippocampus,
    EIR_left_hippo = Left.Hippocampus,
    sex = Gender,
    age = `Age..at.imaging.`,
    APOE = APOE,
    E4_binary = E4,
    CDR = CDR,
    PET_compound = `PET.Compound`,
    PET_fSUVR_rsf_TOT_CORTMEAN = PET_fSUVR_rsf_TOT_CORTMEAN,
    SUVR_threshold = `SUVR.Threshold`,
    amyloid_binary = `Amyloid.status`,
    APOE_npscore = APOE_npscore,
    CDR_binary_point5 = CDR_binary,
    dupes_ordered = `dupes..in.order.`
    # ROI columns "X1" to "X105" remain the same
  )

## Saved merged_df
oasis_df_all_path <- glue("{data_path}/preprocessed/oasis_df_all_apoenpscore_CDRbinary_dupes_allregions_allfuncnetworks_DMNaggregated.csv")
write.csv(df, oasis_df_all_path)




#-###########################-#
## Create the longitudinal df (Similar steps)  --
#-###########################-#

df_longitudinal <- df %>%
  group_by(subject_ID) %>%
  filter(n() >= 3) %>%            # Subjects with at least 3 scans
  filter(any(CDR == 0)) %>%       # Subjects with at least one CDR = 0
  ungroup() %>%
  droplevels()

# Remove rows with missing values for the necessary variables
var_resp <- "EIR_whole_brain"
var_fac <- c("sex", "E4_binary", "APOE", "CDR", "CDR_binary_point5")
var_num <- c("age")
var_rand <- "subject_ID"

df_longitudinal <- df_longitudinal %>%
  drop_na(all_of(c(var_num, var_fac, var_rand, var_resp)))

# Center age by age at first scan
df_longitudinal_agebysubjectcentered <- df_longitudinal %>%
  group_by(subject_ID) %>%
  mutate(
    years_since_first_scan = age - min(age, na.rm = TRUE),
    age_at_first_scan = min(age, na.rm = TRUE)
  ) %>%
  ungroup()

# Adjust PET variables if needed
df_longitudinal_agebysubjectcentered$fSUVR <- df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN
df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN <- NULL




#-###########################-#
## Run LMM for Aggregated DMN and Limbic Func Netowrks EIR Columns ----
#-###########################-#

# Identify all functional network columns
#OLD CODE: all_func_cols <- names(df_longitudinal_agebysubjectcentered)[grepl("^func__", names(df_longitudinal_agebysubjectcentered))]
func_cols <- grep("^func__", names(df_longitudinal_agebysubjectcentered), value = TRUE)
df_longitudinal_agebysubjectcentered[func_cols] <- lapply(df_longitudinal_agebysubjectcentered[func_cols], function(x) as.numeric(as.character(x)))

# EIR types we consider
eir_types <- c("in", "between", "in_between")

all_results <- list()

dir_results_funcnetworks = paste0(results_dir, "/func_networks")

for (eir_type in eir_types) {
  # Identify columns for this EIR type
  type_cols <- func_cols[grepl(paste0("^func__", eir_type, "__"), func_cols)]
  
  # Initialize a data frame to store results for this type
  results_df <- data.frame()
  
  for (col_name in type_cols) {
    # Check if dependent variable contains Inf, -Inf, or NaN
    dep_values <- df_longitudinal_agebysubjectcentered[[col_name]]
    if (any(is.infinite(dep_values) | is.nan(dep_values))) {
      message(glue("Dependent variable is inf, -inf, or NaN for {col_name}, so skipping this test"))
      next
    }
    
    # Extract the network number from the column name
    # Assumes format: func__{eir_type}__{net}
    parts <- strsplit(col_name, "__")[[1]]
    net_num <- parts[length(parts)]
    
    # Build the formula (same structure as before)
    formula_str <- paste0("`", col_name, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")
    model_formula <- as.formula(formula_str)
    
    model <- tryCatch({
      lmer(model_formula, data = df_longitudinal_agebysubjectcentered)
    }, error = function(e) {
      message(glue("Error in fitting model for {col_name}: {e$message}"))
      return(NULL)
    })
    
    if (is.null(model)) {
      next
    }
    
    model_summary <- summary(model)
    coefs <- as.data.frame(model_summary$coefficients)
    
    # Add Network and EIR type info
    coefs$EIR_type <- eir_type
    coefs$Network <- net_num
    coefs$Term <- rownames(coefs)
    
    # Rename columns to match style
    colnames(coefs)[colnames(coefs) == "Std. Error"] <- "Std_Error"
    colnames(coefs)[colnames(coefs) == "t value"] <- "t_value"
    colnames(coefs)[colnames(coefs) == "Pr(>|t|)"] <- "p_value"
    
    coefs <- coefs[, c("EIR_type", "Network", "Term", "Estimate", "Std_Error", "df", "t_value", "p_value")]
    results_df <- rbind(results_df, coefs)
  }
  
  # Remove rows with NA p-values
  results_df <- results_df[!is.na(results_df$p_value), ]
  
  # Apply FDR correction
  results_df$FDR_p_value <- p.adjust(results_df$p_value, method = "fdr")
  
  # Save results for this EIR type
  all_results[[eir_type]] <- results_df
  
  # You can write out CSVs if desired:
  write.csv(results_df, glue("{dir_results_funcnetworks}/functional_networks_DMNaggregated_all_terms_{eir_type}.csv"), row.names = FALSE)
  
  # Significant raw p-values
  signif_results <- results_df %>% filter(p_value < 0.05)
  write.csv(signif_results, glue("{dir_results_funcnetworks}/functional_networks_DMNaggregated_signif_raw_{eir_type}.csv"), row.names = FALSE)
  
  # Significant FDR
  signif_results_FDR <- results_df %>% filter(FDR_p_value < 0.05)
  write.csv(signif_results_FDR, glue("{dir_results_funcnetworks}/functional_networks_DMNaggregated_signif_FDR_{eir_type}.csv"), row.names = FALSE)
}

# Combine all types into one final results dataframe if needed
final_results <- do.call(rbind, all_results)
write.csv(final_results, glue("{dir_results_funcnetworks}/functional_networks_DMNaggregated_all_EIR_types_results.csv"), row.names = FALSE)




#-###########################-#
## Plot results just for DMN (in) ----
#-###########################-#

# Define column name
col_name <- "func__in__23"

# Build the formula (same structure as before)
formula_str <- paste0("`", col_name, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")
model_formula <- as.formula(formula_str)

model <- tryCatch({
  lmer(model_formula, data = df_longitudinal_agebysubjectcentered)
}, error = function(e) {
  message(glue("Error in fitting model for {col_name}: {e$message}"))
  return(NULL)
})

# Obtain estimated marginal means for your binary variable
em <- emmeans(model, specs = ~ E4_binary)

# Convert emmeans object to a data frame for plotting
em_df <- as.data.frame(em)

# em_df will contain:
# - binary_var: the levels of your binary variable
# - emmean: the predicted mean of the outcome at that level
# - SE, df, lower.CL, upper.CL: standard error, degrees of freedom, and confidence limits

# Get the contrast (difference) between the levels of binary_var
contrast_df <- summary(contrast(em, "pairwise"))

# Plot the predicted means with 95% CI error bars
ggplot(em_df, aes(x = E4_binary, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(title = "Predicted Means by Binary Variable",
       x = "Binary Variable",
       y = "Predicted Mean Outcome") +
  theme_minimal()

# Plot with customizations
custom_colors <- c("orange", "black")
plot <- ggplot(em_df, aes(x = E4_binary, y = emmean, color = E4_binary)) +
               geom_point(size = 5) +  # Enlarged points
               geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 1) +  # Error bars
               geom_text(aes(label = round(emmean, 2)), vjust = -1, hjust = 1.2, size = 5) +  # Mean values displayed above points
               scale_color_manual(values = custom_colors) +
               theme_minimal() +
               theme(
                 axis.line.x = element_line(size = 0.5, color = "black"),  # Add bottom axis line
                 axis.line.y = element_line(size = 0.5, color = "black"),  # Add left axis line
                 axis.ticks.x = element_blank(),  # Remove x-axis tick marks
                 axis.text.x = element_blank(),  # Remove x-axis tick labels
                 axis.text.y = element_text(size=16),
                 axis.title.x = element_blank(),  # Remove x-axis title
                 axis.title.y = element_blank(),  # Remove y-axis title
                 panel.grid.major.y = element_line(color = "grey90", size = 0.2),  # Add faint y-axis gridlines
                 panel.grid.minor.y = element_blank(),
                 panel.grid.major.x = element_blank(),  # Remove x-axis gridlines
                 legend.position = "none",  # Hide legend
                 plot.title = element_text(size = 18, face = "bold")
               )


# Save plot
filepath <- "/Users/aburns/Codebook/Projects/tda_criticality/images/OASIS_Longitudinal_Images/manuscript_images_FINAL/DMN_FuncNetowrk_E4_emmeansplot.pdf"
ggsave(filename = filepath,
       plot = plot,
       width = 4.5,
       height = 6,
       dpi = 300,
       units = "in")





#-#############################################################################################################################################-#




#-###########################-#
# DMN & Limbic: Import df and Preprocess  ----
#-###########################-#

## Import data frame 
oasis_df_path <- glue("{data_path}/preprocessed/oasis_df_all_apoenpscore_CDRbinary_dupes_allregions.csv")
df <- read.csv(oasis_df_path, header=T)[,-1]

## import other dataframes with functional info to merge with above df
filepath_func_df <- paste0(data_path, "/OASIS_Data_Func_Networks_DMNandLimbicaggregated_EIR.xlsx")
func_df = read.xlsx(filepath_func_df)


### Merge func_df with demographic df
func_df$subject_idx = func_df$subject_idx + 1
colnames(df)[colnames(df) == "SUBJECT."] <- "subject_idx"
df <- merge(df, func_df, by="subject_idx", all.x = TRUE)

## Rename columns
df <- df %>%
  dplyr::rename(
    MR_ID = `MR.ID`,
    subject_ID = Subject,
    EIR_whole_brain = Whole.Brain,
    EIR_right_hippo = Right.Hippocampus,
    EIR_left_hippo = Left.Hippocampus,
    sex = Gender,
    age = `Age..at.imaging.`,
    APOE = APOE,
    E4_binary = E4,
    CDR = CDR,
    PET_compound = `PET.Compound`,
    PET_fSUVR_rsf_TOT_CORTMEAN = PET_fSUVR_rsf_TOT_CORTMEAN,
    SUVR_threshold = `SUVR.Threshold`,
    amyloid_binary = `Amyloid.status`,
    APOE_npscore = APOE_npscore,
    CDR_binary_point5 = CDR_binary,
    dupes_ordered = `dupes..in.order.`
    # ROI columns "X1" to "X105" remain the same
  )

## Saved merged_df
oasis_df_all_path <- glue("{data_path}/preprocessed/oasis_df_all_apoenpscore_CDRbinary_dupes_allregions_allfuncnetworks_DMNaggregated.csv")
write.csv(df, oasis_df_all_path)





#-###########################-#
## Create the longitudinal df (Similar steps)  ----
#-###########################-#

df_longitudinal <- df %>%
  group_by(subject_ID) %>%
  filter(n() >= 3) %>%            # Subjects with at least 3 scans
  filter(any(CDR == 0)) %>%       # Subjects with at least one CDR = 0
  ungroup() %>%
  droplevels()

# Remove rows with missing values for the necessary variables
var_resp <- "EIR_whole_brain"
var_fac <- c("sex", "E4_binary", "APOE", "CDR", "CDR_binary_point5")
var_num <- c("age")
var_rand <- "subject_ID"

df_longitudinal <- df_longitudinal %>%
  drop_na(all_of(c(var_num, var_fac, var_rand, var_resp)))

# Center age by age at first scan
df_longitudinal_agebysubjectcentered <- df_longitudinal %>%
  group_by(subject_ID) %>%
  mutate(
    years_since_first_scan = age - min(age, na.rm = TRUE),
    age_at_first_scan = min(age, na.rm = TRUE)
  ) %>%
  ungroup()

# Adjust PET variables if needed
df_longitudinal_agebysubjectcentered$fSUVR <- df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN
df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN <- NULL






#-###########################-#
## Run LMM for Aggregated DMN and Limbic Func Netowrks EIR Columns ----
#-###########################-#

# Identify all functional network columns
#OLD CODE: all_func_cols <- names(df_longitudinal_agebysubjectcentered)[grepl("^func__", names(df_longitudinal_agebysubjectcentered))]
func_cols <- grep("^func__", names(df_longitudinal_agebysubjectcentered), value = TRUE)
df_longitudinal_agebysubjectcentered[func_cols] <- lapply(df_longitudinal_agebysubjectcentered[func_cols], function(x) as.numeric(as.character(x)))

# EIR types we consider
eir_types <- c("in", "between", "in_between")

all_results <- list()

dir_results_funcnetworks = paste0(results_dir, "/func_networks")

for (eir_type in eir_types) {
  # Identify columns for this EIR type
  type_cols <- func_cols[grepl(paste0("^func__", eir_type, "__"), func_cols)]
  
  # Initialize a data frame to store results for this type
  results_df <- data.frame()
  
  for (col_name in type_cols) {
    # Check if dependent variable contains Inf, -Inf, or NaN
    dep_values <- df_longitudinal_agebysubjectcentered[[col_name]]
    if (any(is.infinite(dep_values) | is.nan(dep_values))) {
      message(glue("Dependent variable is inf, -inf, or NaN for {col_name}, so skipping this test"))
      next
    }
    
    # Extract the network number from the column name
    # Assumes format: func__{eir_type}__{net}
    parts <- strsplit(col_name, "__")[[1]]
    net_num <- parts[length(parts)]
    
    # Build the formula (same structure as before)
    formula_str <- paste0("`", col_name, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")
    model_formula <- as.formula(formula_str)
    
    model <- tryCatch({
      lmer(model_formula, data = df_longitudinal_agebysubjectcentered)
    }, error = function(e) {
      message(glue("Error in fitting model for {col_name}: {e$message}"))
      return(NULL)
    })
    
    if (is.null(model)) {
      next
    }
    
    model_summary <- summary(model)
    coefs <- as.data.frame(model_summary$coefficients)
    
    # Add Network and EIR type info
    coefs$EIR_type <- eir_type
    coefs$Network <- net_num
    coefs$Term <- rownames(coefs)
    
    # Rename columns to match style
    colnames(coefs)[colnames(coefs) == "Std. Error"] <- "Std_Error"
    colnames(coefs)[colnames(coefs) == "t value"] <- "t_value"
    colnames(coefs)[colnames(coefs) == "Pr(>|t|)"] <- "p_value"
    
    coefs <- coefs[, c("EIR_type", "Network", "Term", "Estimate", "Std_Error", "df", "t_value", "p_value")]
    results_df <- rbind(results_df, coefs)
  }
  
  # Remove rows with NA p-values
  results_df <- results_df[!is.na(results_df$p_value), ]
  
  # Apply FDR correction
  results_df$FDR_p_value <- p.adjust(results_df$p_value, method = "fdr")
  
  # Save results for this EIR type
  all_results[[eir_type]] <- results_df
  
  # You can write out CSVs if desired:
  write.csv(results_df, glue("{dir_results_funcnetworks}/functional_networks_DMNandLIMBICaggregated_all_terms_{eir_type}.csv"), row.names = FALSE)
  
  # Significant raw p-values
  signif_results <- results_df %>% filter(p_value < 0.05)
  write.csv(signif_results, glue("{dir_results_funcnetworks}/functional_networks_DMNandLIMBICaggregated_signif_raw_{eir_type}.csv"), row.names = FALSE)
  
  # Significant FDR
  signif_results_FDR <- results_df %>% filter(FDR_p_value < 0.05)
  write.csv(signif_results_FDR, glue("{dir_results_funcnetworks}/functional_networks_DMNandLIMBICaggregated_signif_FDR_{eir_type}.csv"), row.names = FALSE)
}

# Combine all types into one final results dataframe if needed
final_results <- do.call(rbind, all_results)
write.csv(final_results, glue("{dir_results_funcnetworks}/functional_networks_DMNandLIMBICaggregated_all_EIR_types_results.csv"), row.names = FALSE)






#-###########################-#
## Plot results DMN and Limbic (in-between) ----
#-###########################-#

# Define column name
col_name <- "func__in_between__23"

# Build the formula (same structure as before)
formula_str <- paste0("`", col_name, "` ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID)")
model_formula <- as.formula(formula_str)

model <- tryCatch({
  lmer(model_formula, data = df_longitudinal_agebysubjectcentered)
}, error = function(e) {
  message(glue("Error in fitting model for {col_name}: {e$message}"))
  return(NULL)
})

#Obtain estimated marginal means for your binary variable
em <- emmeans(model, specs = ~ sex)

# Convert emmeans object to a data frame for plotting
em_df <- as.data.frame(em)

# em_df will contain:
# - binary_var: the levels of your binary variable
# - emmean: the predicted mean of the outcome at that level
# - SE, df, lower.CL, upper.CL: standard error, degrees of freedom, and confidence limits

# Get the contrast (difference) between the levels of binary_var
contrast_df <- summary(contrast(em, "pairwise"))

# Plot the predicted means with 95% CI error bars
ggplot(em_df, aes(x = sex, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(title = "Predicted Means by Binary Variable",
       x = "Binary Variable",
       y = "Predicted Mean Outcome") +
  theme_minimal()

custom_colors <- c("orange", "black")

# Plot with customizations
plot <- ggplot(em_df, aes(x = sex, y = emmean, color = sex)) +
  
  geom_jitter(data = df_longitudinal_agebysubjectcentered, 
              aes(x = factor(sex), y = !!sym(col_name)), 
              width = 0.05, alpha = 0.8, color = "grey80", size = 1) +
  
  geom_point(size = 5) +  # Enlarged points
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 1) +  # Error bars
  geom_text(aes(label = round(emmean, 3)), vjust = -1, hjust = 1.2, size = 8) +  # Mean values displayed above points
  scale_color_manual(values = custom_colors) +
  #geom_jitter(width = 0.2, alpha = 0.6) + 
  theme_minimal() +
  theme(
    axis.line.x = element_line(size = 0.5, color = "black"),  # Add bottom axis line
    axis.line.y = element_line(size = 0.5, color = "black"),  # Add left axis line
    axis.ticks.x = element_blank(),  # Remove x-axis tick marks
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_text(size=24, color = "black"),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    panel.grid.major.y = element_line(color = "grey80", size = 0.2),  # Add faint y-axis gridlines
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),  # Remove x-axis gridlines
    legend.position = "none",  # Hide legend
    plot.title = element_text(size = 18, face = "bold")
  )

print(plot)

# Save plot
filepath <- "/Users/aburns/Codebook/Projects/tda_criticality/images/OASIS_Longitudinal_Images/manuscript_images_FINAL/DMNandLimbic_FuncNetowrk_sex_emmeansplot_jitteredpoints.pdf"
ggsave(filename = filepath,
       plot = plot,
       width = 6,
       height = 4.5,
       dpi = 300,
       units = "in")




#### Decided not to plot below because it is expected.

### Now let's do it for years_since_fist_scan
# For convenience, rename the data
df2 <- df_longitudinal_agebysubjectcentered

# Calculate 5th and 95th percentile
age_low  <- quantile(df2$years_since_first_scan, probs = 0.05, na.rm = TRUE)
age_high <- quantile(df2$years_since_first_scan, probs = 0.95, na.rm = TRUE)

# Create a sequence of e.g. 50 values between the 5th and 95th percentile
age_seq <- seq(age_low, age_high, length.out = 50)

# emtrends() calculates how the outcome changes per unit change
# in the variable you specify (here, "scale(age_at_first_scan)"),
# separately for each level of E4_binary.
age_trends <- emtrends(
  model,          
  var = "years_since_first_scan"  # the slope you want
)

# Compare the difference in slopes between E4 = 0 vs. E4 = 1
contrast_results <- contrast(age_trends, method = "pairwise")

summary(contrast_results)

# emmeans for the interaction of age_at_first_scan by E4_binary
# "type" can be "response" if you want predictions on the outcome scale (for linear models it's the same as "link")
em_interaction <- emmeans(
  model,
  specs = ~ age_at_first_scan | E4_binary,
  at = list(age_at_first_scan = age_seq),
  type = "response"
)

# Convert emmeans results to a data frame
df_interaction <- as.data.frame(em_interaction)








#-#######################-#
# Export Source data for manuscript plot ----
#-#######################-#
# Figure 4 is made from df_longitudinal_agebysubjectcentered

#Filter df to make sure it abides by HIPPA
df_fig4sourcedata <- df_longitudinal_agebysubjectcentered %>%
  dplyr::select(-c("subject_idx", "dupes_ordered")) %>%
  mutate(across(-MR_ID, ~ subject_ID))

# Export data to csv
source_data_path <- "/Users/aburns/Codebook/Projects/tda_criticality/papers_and_protocols/MANUSCRIPT_OASIS_Longitudinal/CommBio_Submission/Accepted Version/SourceData_Figure4_DMNLimbic.csv"
write.csv(df_longitudinal_agebysubjectcentered, source_data_path)
source_data_path <- "/Users/aburns/Codebook/Projects/tda_criticality/papers_and_protocols/MANUSCRIPT_OASIS_Longitudinal/CommBio_Submission/Accepted Version/SourceData_Figure4_DMNLimbic_OASIScompliant.csv"
write.csv(df_fig4sourcedata, source_data_path)

















### OLD Code for age at first scan and E4 for just "in"

# For convenience, rename the data
df2 <- df_longitudinal_agebysubjectcentered

# Calculate 5th and 95th percentile
age_low  <- quantile(df2$age_at_first_scan, probs = 0.05, na.rm = TRUE)
age_high <- quantile(df2$age_at_first_scan, probs = 0.95, na.rm = TRUE)

# Create a sequence of e.g. 50 values between the 5th and 95th percentile
age_seq <- seq(age_low, age_high, length.out = 50)

# emtrends() calculates how the outcome changes per unit change
# in the variable you specify (here, "scale(age_at_first_scan)"),
# separately for each level of E4_binary.
age_trends <- emtrends(
  model,
  ~ E4_binary,            # grouping factor whose levels you'll compare
  var = "age_at_first_scan"  # the slope you want
)

# Compare the difference in slopes between E4 = 0 vs. E4 = 1
contrast_results <- contrast(age_trends, method = "pairwise")

summary(contrast_results)

# emmeans for the interaction of age_at_first_scan by E4_binary
# "type" can be "response" if you want predictions on the outcome scale (for linear models it's the same as "link")
em_interaction <- emmeans(
  model,
  specs = ~ age_at_first_scan | E4_binary,
  at = list(age_at_first_scan = age_seq),
  type = "response"
)

# Convert emmeans results to a data frame
df_interaction <- as.data.frame(em_interaction)

custom_colors <- c("green", "violet")  # for example

plot_interaction <- ggplot(
  df_interaction,
  aes(x = age_at_first_scan, y = emmean, color = E4_binary)
) +
  # 1) Lines for predicted means
  geom_line(size = 1.2) +
  
  # 2) Ribbons for confidence intervals
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL, fill = E4_binary),
    alpha = 0.2,     # opacity
    color = NA       # no outline
  ) +
  
  # 3) Manual color/fill scales, if desired
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # 4) Title, axis labels, and theme
  labs(
    title = "Interaction of Age at First Scan & E4 on Outcome",
    x = "Age at First Scan",
    y = "Predicted Outcome"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

# Print or save
print(plot_interaction)

filepath <- "/Users/aburns/Codebook/Projects/tda_criticality/images/OASIS_Longitudinal_Images/manuscript_images_FINAL/DMN_FuncNetowrk_Age_by_E4_interaction.pdf"
ggsave(
  filepath,
  plot = plot_interaction,
  width = 5,
  height = 4,
  dpi = 300
)












