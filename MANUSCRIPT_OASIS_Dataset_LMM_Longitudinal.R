#---------------------------------------------------#
#  The following code was creted by Andrew P. Burns
#  
#  Note: Much of the workflow followed the template of 
#  the glmmTMB package and pipeline. Please consult 
#  their paper for further clarification.
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
font_add("Arial", regular= "arial.ttf")
showtext_auto()


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








#-###########################-#
# Import df and Preprocess  ----
#-###########################-#

## Import data frame (example only. Any other import function of R can be used).
oasis_df_path <- glue("{data_path}/preprocessed/oasis_df_all_apoenpscore_CDRbinary_dupes.csv")
df <- read.csv(oasis_df_path, header=T)[,-1]


## Carefully check data structure, column names and vector classes. Change them as needed.
str(df)
factor_columns <- sapply(df, is.factor)


##Convert column names
col_names <- c("data_instance", "subject_ID", "MR_ID", "EIR_whole_brain", "EIR_right_hippo", "EIR_left_hippo", "sex",
               "age", "APOE","E4_binary", "CDR", "PET_compound", "PET_fSUVR_rsf_TOT_CORTMEAN", "SUVR_threshold", "amyloid_binary",
               "APOE_npscore", "CDR_binary_point5", "dupes_ordered")
colnames(df) <- col_names


## Add CDR binary for 1
df <- df %>%
  mutate(CDR_binary_1 = if_else(CDR >= 1, 1, 0))


## Convert categorical columns to factors
col_factors <- c("data_instance", "subject_ID", "MR_ID", "sex", "APOE", "E4_binary", "amyloid_binary", "CDR_binary_point5", "CDR_binary_1","PET_compound")
for (col in col_factors){
  df[[col]] <- factor(df[[col]])
}

## Convert numerical CDR to an ordered factor
df$CDR <- factor(df$CDR, ordered=T, levels=c(0,0.5,1,2))


## Add Hyperexcitation Index (# positive / total) for whole_brain, right_hippo, and left_hippo
## For 132 regions, there are 8646 connections (whole_brain). For 105 regions there are 5460 connections (whole_brain)
## For each region, there are either 131 or 104 connections
## Formula for HI from EIR is "(num_connections * EIR) / (1 + EIR)"
num_region_connections <- 104
num_wholebrain_connections <- 8646
df$HI_whole_brain <- (df$EIR_whole_brain * num_wholebrain_connections) / (1 + df$EIR_whole_brain)


## Make longitudinal df with just subjects with >= 3 scans
df_longitudinal <- df %>%
  group_by(subject_ID) %>%
  filter(max(dupes_ordered) >= 3) %>%
  ungroup() %>%
  droplevels()


## Make cross-sectional and longitudinal df with only PET patients
df_PET <- df %>%
  filter(PET_compound != "No PET")

df_longitudinal_PET <- df_PET %>%
  group_by(subject_ID) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  droplevels()









#-##############################-#
# Definition of variables  ----
#-##############################-#

#================================#
# * ONE RESPONSE variable ----
#================================#
var_resp <- "EIR_whole_brain"


#=============================================================================#
# *  Fixed predictors: quantitative and categorical predictor variables ----
#=============================================================================#
## FACTOR PREDICTOR variable(s)
var_fac <- c("sex", "E4_binary", "APOE", "CDR", "CDR_binary_point5","CDR_binary_1")
var_fac_PET <- c("sex", "E4_binary", "APOE", "CDR", "CDR_binary_point5", "CDR_binary_1", "amyloid_binary")


## NUMERIC or INTEGER PREDICTOR variable(s) 
var_num <- c("age")                              
var_num_PET <- c("PET_fSUVR_rsf_TOT_CORTMEAN", "age")                              


#==============================================================#
# * Random predictors: dependency structure of the data ----
#==============================================================#
## RANDOM term(s)
var_rand <- "subject_ID"   


#==============================================#
# * Temporal and spatial data structure ----
#==============================================#
# The variables specified here will be used to check for temporal and/or spatial autocorrelation in model residuals (6.3).

# NUMERIC or INTEGER TIME variable that specifies temporal structure.
var_time <- "age"                          # assign default NA if missing.

# Second, if these temporal data are structured, add the grouping variable below.
var_time_groups <- "subject_ID"                  # assign NA if missing.

# NUMERIC or INTEGER COORDINATES (x, y) that specify spatial structure.
var_space <- NA                         # assign NA if missing.



#=========================#
# * Missing values ----
#=========================#
# We prune the dataset to complete cases for all variables that are considered for analysis.
# NOTE-1: Data pruning should be kept to the minimum needed to allow models to run. 
#         Therefore, rerun 1. after a final model formulation has been identified. This may allow to rescue some observations.
# NOTE-2: A key assumption is that missing data (and the reasons for data to be missing) are randomly distributed across the dataset. 
#         Confirm by inspecting the removed data rows after creating the dataset df.NAs 
df_longitudinal.pr <- remove_NAs(data = df_longitudinal, variables = c(var_num, var_fac, var_rand, var_resp))
df_longitudinal.NAs <- anti_join(df_longitudinal,  df_longitudinal.pr) # This object contains all removed data rows.
df_longitudinal <- df_longitudinal.pr   # Keep only complete cases in the dataset for further analysis: 

# Repeat for PET longitudinal dataframe
df_longitudinal_PET.pr <- remove_NAs(data = df_longitudinal_PET, variables = c(var_num_PET, var_fac_PET, var_rand, var_resp))
df_longitudinal_PET.NAs <- anti_join(df_longitudinal_PET,  df_longitudinal_PET.pr) # This object contains all removed data rows.
df_longitudinal_PET <- df_longitudinal_PET.pr








#-###############################-#
#  Raw data exploration ----
#-###############################-#
# This step helps to identify patterns that may require consideration for model formulation.
# NOTE: Further checks can be found in the residual plot analysis (section 6.1), in particular for missing predictors or interactions.
## Create directory to save plots and results
raw_data_exploration_dir <- glue("{results_dir}/raw_data_exploration")
if (!file.exists(raw_data_exploration_dir)) {dir.create(raw_data_exploration_dir, recursive = TRUE)}


#=========================#
# * Extreme values ----
#=========================#

# We graphically inspect variables for extreme values.
#---------------------------------------------#
# **  Extremes in NUMERICAL variables ---- 
#---------------------------------------------#
## Create directory to save plots and results
extreme_val_dir <- glue("{raw_data_exploration_dir}/extreme_vals_NUMERICAL")
if (!file.exists(extreme_val_dir)) {dir.create(extreme_val_dir, recursive = TRUE)}

## Dotplot of all NUMERICAL variables (predictors + response) WITH and WITHOUT PET
#Whole Brain
plot <- dotplot_num(data = df_longitudinal, variables = c(var_num, var_resp)) + 
  theme(strip.text.x = element_text(size = 12)) +
  labs(title="OASIS Longitudinal (>=3) All")
ggsave(glue("{extreme_val_dir}/dotplot_NUMERICAL_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

plot <- dotplot_num(data = df_longitudinal_PET, variables = c(var_num_PET, var_resp)) + 
  theme(strip.text.x = element_text(size = 8)) +
  labs(title="OASIS Longitudinal (>=3) With PET")
ggsave(glue("{extreme_val_dir}/dotplot_NUMERICAL_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

# What should I look for?
# >> Is any observation CLEARLY separated from the core distribution? Such observersation may represent 'implausible extremes'

# Resolving implausible extreme values: 
# >> Check where such values originate from. If you can trace them to objective (!) typing errors, correct these values.
# >> If extremes are in the response variable: choose an adequate distribution family.
# >> If extremes are in predictor variables: 
#    1. Use model assessment to check if these predictor values cause concern.
#    2. If so, consider data transformation to mitigate the issue.
# >> If extremes cannot be modeled appropriately, consider reporting effect estimates and their SE with and without these extremes.


## Histograms of all NUMERICAL variables (predictors + response)
binwidths <- setNames(c(2, 0.05), c(var_num, var_resp))
plot <- histogram_num_drew(data = df_longitudinal, variables = c(var_num, var_resp), binwidths = binwidths) + 
  theme(strip.text.x = element_text(size = 8)) +
  labs(title="OASIS Longitudinal (>=3) All")
ggsave(glue("{extreme_val_dir}/histograms_NUMERICAL_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

binwidths <- setNames(c(0.1, 2, 0.05), c(var_num_PET, var_resp))
plot <- histogram_num_drew(data = df_longitudinal_PET, variables = c(var_num_PET, var_resp), binwidths=binwidths) + 
  theme(strip.text.x = element_text(size = 8)) +
  labs(title="OASIS Longitudinal (>=3) With PET")
ggsave(glue("{extreme_val_dir}/histograms_NUMERICAL_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)


## Make histogram of PET fSUVR by PET tracer
PIB_thresh <- 1.42
AV45_thresh <- 1.19
plot <- ggplot(df_longitudinal_PET, aes(x = PET_fSUVR_rsf_TOT_CORTMEAN, fill = PET_compound)) +
  geom_histogram(data = subset(df_longitudinal_PET, PET_compound == "PIB"), alpha = 0.5, position = "identity", binwidth = 0.1) +
  geom_histogram(data = subset(df_longitudinal_PET, PET_compound == "AV45"), alpha = 0.5, position = "identity", binwidth = 0.1) +
  geom_vline(xintercept = PIB_thresh, linetype = "dashed", color = "blue", size = 1) +
  geom_vline(xintercept = AV45_thresh, linetype = "dashed", color = "red", size = 1) +
  scale_fill_manual(values = c("PIB" = "blue", "AV45" = "red")) +
  labs(x = "fSUVR", y = "Count", fill = "PET Tracer") +
  theme_minimal()
ggsave(glue("{extreme_val_dir}/histogram_PETcompounds.pdf"), plot = plot, width = 10, height = 6, dpi = 300)



#------------------------------------------#  
  # ** Extremes in FACTOR variables ----      
#------------------------------------------#
## Create directory to save plots and results
extreme_val_dir <- glue("{raw_data_exploration_dir}/extreme_vals_FACTOR")
if (!file.exists(extreme_val_dir)) {dir.create(extreme_val_dir, recursive = TRUE)}

## Plot barplots
plot <- barplot_fac(data = df_longitudinal, variables = c(var_fac, var_rand)) +
  theme(strip.text.x = element_text(size = 8)) +
  labs(title="OASIS Longitudinal (>=3) All")
ggsave(glue("{extreme_val_dir}/barplot_FACTOR_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

plot <- barplot_fac(data = df_longitudinal_PET, variables = c(var_fac_PET, var_rand)) +
  theme(strip.text.x = element_text(size = 8)) +
  labs(title="OASIS Longitudinal (>=3) With PET")
ggsave(glue("{extreme_val_dir}/barplot_FACTOR_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

# What should I look for?
# >> Are observations roughly balanced between grouping levels (conforming to a "balanced design")?

# Resolving extreme unbalance: 
# >> If some factor levels arise from typing errors: Correct them.
# >> If there is extreme imbalance in sample size among factor levels, 
#    accept that effect estimates for those with poor replication will associate with large uncertainty.
#    Alternatively, consider pooling levels with very few replicates, but only if this remains biologically meaningful. 



#------------------------------------------#  
# ** Followup on extreme values ----      
#------------------------------------------#
## Followup on extreme values
# 3 subjects with extreme values in EIB_whole_brain
subject_IDs_xtreme_EIB_whole_brain <- df_longitudinal %>%
  filter(EIR_whole_brain >= 1.5) %>%
  dplyr::select(subject_ID) %>%
  pull()

# 1 subject with extreme values in EIB_left_hippo
subject_IDs_xtreme_EIB_hippoL <- df_longitudinal %>%
  filter(EIR_left_hippo >= 2.75) %>%
  dplyr::select(subject_ID) %>%
  pull()

# 1 subject with extreme values in EIB_right_hippo
subject_IDs_xtreme_EIB_hippoR <- df_longitudinal %>%
  filter(EIR_right_hippo >= 2.45) %>%
  dplyr::select(subject_ID) %>%
  pull()

# 2 subjects with extreme values in SUVR
xtreme_subjs <- reduce(list(subject_IDs_xtreme_EIB_whole_brain, subject_IDs_xtreme_EIB_hippoL, subject_IDs_xtreme_EIB_hippoR), union)
xtreme_subjs_PET <- intersect(xtreme_subjs, df_longitudinal_PET$subject_ID)


xtreme_colors <- c("red", "blue", "green", "orange", "purple")
xtreme_colors_PET <- c("red", "blue", "green")
names(xtreme_colors) <- xtreme_subjs
names(xtreme_colors_PET) <- xtreme_subjs_PET

plot <- dotplot_num_highlight_subjects(data = df_longitudinal, variables = c(var_num, var_resp), highlight_subjects = xtreme_colors) +
  theme(strip.text.x = element_text(size = 10)) +
  labs(title="OASIS Longitudinal (>=3): Extreme Values (colored)")
ggsave(glue("{extreme_val_dir}/dotplot_xtreme_NUMERICAL_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

plot <- dotplot_num_highlight_subjects(data = df_longitudinal_PET, variables = c(var_num_PET, var_resp), highlight_subjects = xtreme_colors_PET) + #Subject ID OAS30062 not in PET df
  theme(strip.text.x = element_text(size = 10)) +
  labs(title="OASIS Longitudinal (>=3) With PET: Extreme Values (colored)")
ggsave(glue("{extreme_val_dir}/dotplot_xtreme_NUMERICAL_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)




#================================#
# * Predictor collinearity ----
#================================#
# We graphically and numerically inspect variables for predictor collinearity.
collinearity_dir <- glue("{raw_data_exploration_dir}/collinearity_and_vif")
if (!file.exists(collinearity_dir)) {dir.create(collinearity_dir, recursive = TRUE)}

#-------------------------------------------------------------#
# **  Graphical inspection for predictor collinearity ----
#-------------------------------------------------------------#
# We graphically inspect predictor collinearity for all pairwise combinations of numeric and/or factor predictors.

# ***  NUMERIC predictors: Scatterplots ----
# Skip this step if you have < 2 numeric predictors
# What should I look for?
# >> Data should distribute rather homogeneously.


# *** FACTOR against NUMERIC predictors: Swarm plots ----
# Skip this step if your predictors are either all numeric, OR all factor.
# Swarm boxplots for all numeric against all factor predictors.
plot <- coll_num_fac_drew(data = df_longitudinal, predictors = c(var_num, var_fac), title_font_size = 8, x_axis_font_size = 8, y_axis_font_size = 8) 
ggsave(glue("{collinearity_dir}/collinearity_NUMERICALvsFACTOR_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

plot <- coll_num_fac_drew(data = df_longitudinal_PET, predictors = c(var_num_PET, var_fac_PET), title_font_size = 10, x_axis_font_size = 8, y_axis_font_size = 8)
ggsave(glue("{collinearity_dir}/collinearity_NUMERICALvsFACTOR_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

# What should I look for?
# >> Numeric values on the y axis should distribute roughly homogeneously across levels on the x-axis.


# *** FACTOR predictors: Mosaic plots ----   
# Skip this step if you have < 2 factor predictors.
plot <- coll_fac_drew(data = df_longitudinal, predictors = var_fac, title_font_size = 8, x_axis_font_size = 8, y_axis_font_size = 8)
ggsave(glue("{collinearity_dir}/collinearity_FACTORvsFACTOR_MosaicPlot_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

plot <- coll_fac_drew(data = df_longitudinal_PET, predictors = var_fac_PET, title_font_size = 8, x_axis_font_size = 8, y_axis_font_size = 8)
ggsave(glue("{collinearity_dir}/collinearity_FACTORvsFACTOR_MosaicPlot_WithPET.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

# What should I look for?
# >> All combinations of predictor levels should have roughly similar sample sizes.



#------------------------------------------------#
# ** Variance Inflation Factors (VIFs) ----
# # Derived from Zuur, Ieno & Elphick (2010).
#------------------------------------------------#
#Run VIF on every model predictor combination (both binary and continous versions of predictor)

# No CDR, APOE4 (binary), amyloid (SUVR)
corvif(data = df_longitudinal, variables = c("age", "E4_binary", "sex"))
corvif(data = df_longitudinal_PET, variables = c("age", "E4_binary", "sex", "PET_fSUVR_rsf_TOT_CORTMEAN"))
# No CDR, APOE4 (binary), amyloid (binary)
corvif(data = df_longitudinal_PET, variables = c("age", "E4_binary", "sex", "amyloid_binary"))

# CDR (categorical), APOE4 (binary), amyloid (SUVR)
corvif(data = df_longitudinal, variables = c("age", "CDR", "E4_binary", "sex"))
corvif(data = df_longitudinal_PET, variables = c("age", "CDR", "E4_binary", "sex", "PET_fSUVR_rsf_TOT_CORTMEAN"))
# CDR (categorical), APOE4 (binary), amyloid (binary)
corvif(data = df_longitudinal_PET, variables = c("age", "CDR", "E4_binary", "sex", "amyloid_binary"))

# CDR (binary - 0.5), APOE4 (binary), amyloid (SUVR)
corvif(data = df_longitudinal, variables = c("age", "CDR_binary_point5", "E4_binary", "sex"))
corvif(data = df_longitudinal_PET, variables = c("age","CDR_binary_point5", "E4_binary", "sex", "PET_fSUVR_rsf_TOT_CORTMEAN"))
# CDR (binary - 0.5), APOE4 (binary), amyloid (binary)
corvif(data = df_longitudinal_PET, variables = c("age", "CDR_binary_point5", "E4_binary", "sex", "amyloid_binary"))

# CDR (binary - 1), APOE4 (binary), amyloid (SUVR)
corvif(data = df_longitudinal, variables = c("age", "CDR_binary_1", "E4_binary", "sex"))
#corvif(data = df_longitudinal_PET, variables = c("age","CDR_binary_1", "E4_binary", "sex", "PET_fSUVR_rsf_TOT_CORTMEAN")) #### ERROR: CDR_binary1 only has one level
# CDR (binary - 1), APOE4 (binary), amyloid (binary)
#corvif(data = df_longitudinal_PET, variables = c("age","CDR_binary_1", "E4_binary", "sex", "amyloid_binary")) #### ERROR: CDR_binary1 only has one level


# What should I look for?
# >> Check column "GVIF": Entirely independent predictors yield GVIF = 1, while larger GVIF values indicate 
#   increasing predictor collinearity (see article section 4.2.2 for details)




#==========================================#
#* Predictor-response relationships ----
#==========================================#
predictor_response_dir <- glue("{raw_data_exploration_dir}/predictor_response")
if (!file.exists(predictor_response_dir)) {dir.create(predictor_response_dir, recursive = TRUE)}

# The function generates plots of the response variable against each of the potential predictor variables.
relat_single_drew(data = df_longitudinal, response = var_resp, predictors = c(var_num, var_fac), save_file_path = glue("{predictor_response_dir}/EIR_whole_brain_vs_predictors_all.pdf"))
relat_single_drew(data = df_longitudinal_PET, response = var_resp, predictors = c(var_num_PET, var_fac_PET), save_file_path = glue("{predictor_response_dir}/EIR_whole_brain_vs_predictors_PET.pdf"))

# NOTE: Up to 9 plots are shown in the plot window... 
#       If > 9 plots are present, they are saved in your working directory in a file 
#       named "Predictor_response_relationships.pdf"

# What should I look for?
# >> Identify relevant relationships with the response.




#------------------------------------------------#
# ** Sphaghetti Plots of Predictor-response relationships ----
#------------------------------------------------#
simple_sphaghetti_plot <- function(df, x, y, group, color, shape=NULL, title=glue("{y} vs. {x}")) {
  # Start the ggplot object
  p <- ggplot(df, aes_string(x = x, y = y, group = group, color = color)) +
    geom_line() +
    theme_minimal() +
    labs(title = title, x = x, y = y)
  # Conditionally add geom_point if 'shape' is not NULL
  if (!is.null(shape)) {
    p <- p + geom_point(aes_string(shape = shape))
  }
  # Return the plot object
  return(p)
}

## EIR vs. age (Group=Subject, color=sex)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="EIR_whole_brain", group="subject_ID", color="sex", shape="CDR")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIR_whole_brain_vs_predictors_by_SubjectSex_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="EIR_whole_brain", group="subject_ID", color="E4_binary", shape="CDR")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIR_whole_brain_vs_predictors_by_SubjectE4_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

#Filter just subjects that have a CDR other than 0
filtered_df <- df_longitudinal %>%
  group_by(subject_ID) %>% 
  filter(any(CDR > 0)) %>%  
  ungroup()

## Just subjects that have a CDR other than 0
plot <- simple_sphaghetti_plot(filtered_df, x="age", y="EIR_whole_brain", group="subject_ID", color="sex", shape="CDR") +
  labs(title="Only Subjects with at least one CDR > 0")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIR_whole_brain_vs_predictors_by_SubjectSex_CDRpositive_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)
plot <- simple_sphaghetti_plot(filtered_df, x="age", y="EIR_whole_brain", group="subject_ID", color="E4_binary", shape="CDR") +
  labs(title="Only Subjects with at least one CDR > 0")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIR_whole_brain_vs_predictors_by_SubjectE4_CDR_positive_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

## CDR vs. age (Group=Subject, color=sex)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="CDR", group="subject_ID", color="sex")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIB_CDR_vs_Age_by_SubjectSex_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="CDR", group="subject_ID", color="E4_binary")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIB_CDR_vs_Age_by_SubjectE4_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

## SUVR vs. age (Group=Subject, color=sex)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="PET_fSUVR_rsf_TOT_CORTMEAN", group="subject_ID", color="sex", shape="CDR")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIB_SUVR_vs_Age_by_SubjectSex_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)
plot <- simple_sphaghetti_plot(df_longitudinal, x="age", y="PET_fSUVR_rsf_TOT_CORTMEAN", group="subject_ID", color="E4_binary", shape="CDR")
ggsave(glue("{predictor_response_dir}/SphaghettiPlotEIB_SUVR_vs_Age_by_SubjectE4_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)




#===============================#
# * Response distribution ----
#===============================#
response_dir <- glue("{raw_data_exploration_dir}/response_distributions")
if (!file.exists(response_dir)) {dir.create(response_dir, recursive = TRUE)}

## Run normality tests (Shapiro-Wilk, Kolomogorov-Smirnov, Anderson-Darling, Lilliefors, D-Agostino's K^2 Test)
shapiro_whole_brain <- shapiro.test(df_longitudinal$EIR_whole_brain)
ad_whole_brain <- ad.test(df_longitudinal$EIR_whole_brain)

## QQ plot
plot <- ggqqplot(df_longitudinal$EIR_whole_brain) + labs(title="EIR_whole_brain")
ggsave(glue("{response_dir}/QQPlot_EIR_whole_brain_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)








#-#######################-#
# Model formulation ----
#-#######################-#

#===============================================#
# * Standardisation of numeric predictors ----
#===============================================#

# df: center age by age_at_first_scan
df_agebysubjectcentered <- df
df_agebysubjectcentered <- df_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  mutate(years_since_first_scan = age - min(age, na.rm = TRUE), min_age = min(age, na.rm= T)) %>%
  ungroup()
df$fSUVR <- df_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN 
df_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN <- NULL


# df_longitudinal: center age by age_at_first_scan
df_longitudinal_agebysubjectcentered <- df_longitudinal
df_longitudinal_agebysubjectcentered <- df_longitudinal_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  mutate(years_since_first_scan = age - min(age, na.rm = TRUE), age_at_first_scan = min(age, na.rm= T)) %>%
  ungroup()
df_longitudinal_agebysubjectcentered$fSUVR <- df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN 
df_longitudinal_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN <- NULL


# df_longitudinal_PET: center age by age_at_first_scan
df_longitudinal_PET_agebysubjectcentered <- df_longitudinal_PET
df_longitudinal_PET_agebysubjectcentered <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  mutate(years_since_first_scan = age - min(age, na.rm = TRUE), age_at_first_scan = min(age, na.rm= T)) %>%
  ungroup()
df_longitudinal_PET_agebysubjectcentered$fSUVR <- df_longitudinal_PET_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN 
df_longitudinal_PET_agebysubjectcentered$PET_fSUVR_rsf_TOT_CORTMEAN <- NULL


# Create amyloid negative df. Filter out amyloid for longitudinal PET. Just amyloid - for all >= 3 scans
df_longitudinal_amyloidnegative_agebysubjectcentered <- df_longitudinal_PET_agebysubjectcentered
df_longitudinal_amyloidnegative_agebysubjectcentered <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  filter(amyloid_binary != "Amyloid +") %>% 
  filter(amyloid_binary != "No PET")  %>% 
  group_by(subject_ID) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  droplevels()




#===============================#
# * Model implementation ----
#===============================#


#==============================#
# ** Whole Brain Analysis ----
#==============================#

### Run model with random slope for verification that random slopes don't provide any value
wholebrain_lme4_E4_binary_random_slope <- lmer(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 + scale(years_since_first_scan) | subject_ID), data = df_longitudinal_agebysubjectcentered)
summary_random_slope <- summary(wholebrain_lme4_E4_binary_random_slope)
anova_random_slope <- anova(wholebrain_lme4_E4_binary_random_slope)

### Run model with all longitudinal data using EIR. NO RANDOM SLOPES
wholebrain_lme4_E4_binary <- lmer(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID), data = df_longitudinal_agebysubjectcentered)
summary_binary <- summary(wholebrain_lme4_E4_binary)
anova_binary <- anova(wholebrain_lme4_E4_binary)

### Run model with all longitudinal data using HI
wholebrain_lme4_E4_binary_HI <- lmer(HI_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID), data = df_longitudinal_agebysubjectcentered)
summary_binary_HI <- summary(wholebrain_lme4_E4_binary_HI)
anova_binary_HI <- anova(wholebrain_lme4_E4_binary_HI)

### Run model with amyloid NEGATIVE data 
wholebrain_lme4_E4_binary_amyloidnegative <- lmer(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID), data = df_longitudinal_amyloidnegative_agebysubjectcentered)
summary_amyloidnegative <- summary(wholebrain_lme4_E4_binary_amyloidnegative)
anova_amyloidnegative <- anova(wholebrain_lme4_E4_binary_amyloidnegative)

### Run model with PET data and include amyloid as another covariate: NO INTERACTION
wholebrain_lme4_E4_binary_PET <- lmer(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + amyloid_binary + (1 | subject_ID), data = df_longitudinal_PET_agebysubjectcentered)
summary_PET <- summary(wholebrain_lme4_E4_binary_PET)
anova_PET <- anova(wholebrain_lme4_E4_binary_PET)

### Run model with PET data and include amyloid as another covariate: WITH INTERACTION
wholebrain_lme4_E4_binary_PETinteraction <- lmer(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary * amyloid_binary + (1 | subject_ID), data = df_longitudinal_PET_agebysubjectcentered)
summary_PETinteraction <- summary(wholebrain_lme4_E4_binary_PETinteraction)
anova_PETinteraction <- anova(wholebrain_lme4_E4_binary_PETinteraction)




#==============================#
# ** Save Model summaries and ANOVA summaries ----
#==============================#
sink(file.path(glue("{results_dir}/all_model_anova_results_wholebrain.txt")))

cat("## Random Slope Model Summary\n\n")
print(summary_random_slope)
cat("\n\n## Random Slope Model ANOVA\n\n")
print(anova_random_slope)
cat("\n\n\n\n\n\n\n\n\n\n")

cat("## Binary Model Summary\n\n")
print(summary_binary)
cat("\n\n## Binary Model ANOVA\n\n")
print(anova_binary)
cat("\n\n\n\n\n\n\n\n\n\n")

cat("## Binary HI Model Summary\n\n")
print(summary_binary_HI)
cat("\n\n## Binary HI Model ANOVA\n\n")
print(anova_binary_HI)
cat("\n\n\n\n\n\n\n\n\n\n")

cat("## Amyloid Negative Model Summary\n\n")
print(summary_amyloidnegative)
cat("\n\n## Amyloid Negative Model ANOVA\n\n")
print(anova_amyloidnegative)
cat("\n\n\n\n\n\n\n\n\n\n")

cat("## PET Model Summary\n\n")
print(summary_PET)
cat("\n\n## PET Model ANOVA\n\n")
print(anova_PET)
cat("\n\n\n\n\n\n\n\n\n\n")

cat("## PET Interaction Model Summary\n\n")
print(summary_PETinteraction)
cat("\n\n## PET Interaction Model ANOVA\n\n")
print(anova_PETinteraction)
cat("\n\n\n\n\n\n\n\n\n\n")

sink()




#==============================#
# *** Post-hoc Analysis ----
#==============================#

# Here we decompose the interaction to better understand the slopes
simple_slopes_age_by_sex <- emmeans(wholebrain_lme4_E4_binary, pairwise ~ years_since_first_scan | sex)
simple_slopes_age_by_e4 <- emmeans(wholebrain_lme4_E4_binary, pairwise ~ years_since_first_scan | E4_binary)
simple_slopes_sex_by_e4 <- emmeans(wholebrain_lme4_E4_binary, pairwise ~ sex | E4_binary)
summary(simple_slopes_age_by_sex)
summary(simple_slopes_age_by_e4)
summary(simple_slopes_sex_by_e4)


# Post-hoc comparisons
post_hoc_sex <- emmeans(wholebrain_lme4_E4_binary, pairwise ~ sex | years_since_first_scan * E4_binary)
summary(post_hoc_sex)
post_hoc_e4 <- emmeans(wholebrain_lme4_E4_binary, pairwise ~ E4_binary | years_since_first_scan * sex)
summary(post_hoc_e4)

# Marginal means
marginal_means <- emmeans(wholebrain_lme4_E4_binary, ~ years_since_first_scan * sex * E4_binary)
summary(marginal_means)

# Calculate slopes for age_centered_by_subj across sex and E4_binary
slopes <- emtrends(wholebrain_lme4_E4_binary, pairwise ~ sex * E4_binary, var = "years_since_first_scan")
summary(slopes)

# Extract estimated slopes
estimated_slopes <- emmeans(wholebrain_lme4_E4_binary, ~ sex * E4_binary | years_since_first_scan)
summary(estimated_slopes)






#-#######################-#
# Plotting Results ----
#-#######################-#

#===============================#
# * Set plot settings ----
#===============================#
### Define plotting colors
male_color <- "#c6cb8a"
male_light <- "#f0f3d1"
male_dark <- "#878c45"

female_color <- "#61304f"
female_light <- "#a57a95"
female_dark <- "#2e1826"

age1_color <- "#6A0D0D"
age2_color <- "#F6E8C3"
age3_color <- "#083D4B"

#Color for non-gradient hist
hist_outline_color <- "#011627"
hist_fill_color <- "#A18276"

# Colors for the hist gradient
hist_low_color <- "#B8C6E4"
hist_mid_color <- "dimgrey"
hist_high_color <- "#EB4A26"




#===============================#
# * Plot 3-way interaction term ----
#===============================#

### Replace E4 with epsilon for all df
df_longitudinal_agebysubjectcentered$E4_binary <- factor(df_longitudinal_agebysubjectcentered$E4_binary,
                                                         levels = c("E4-", "E4+"),
                                                         labels = c("\u03B54-", "\u03B54+"))

# Generate the plot for E4- without legend
plot1 <- plot_model(wholebrain_lme4_E4_binary, type = "pred", terms = c("years_since_first_scan", "sex", "E4_binary[\u03B54-]"), axis.labels = "") 
color_levels <- levels(plot1$data$group_col)
color_mapping <- setNames(c(female_color, male_color), color_levels)
plot1 <- plot1 +
  labs(
    title = NULL,
    x = "Years Since First Scan",
    y = "Whole Brain EI Ratio "
  ) +
  scale_color_manual(
    name = "Sex",  # Change this to your desired legend title
    values = color_mapping  # Use the dynamically generated color mapping
  ) +
  scale_fill_manual(
    values = color_mapping  # Use the same color mapping for the confidence intervals
  ) +
  manuscript_theme +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.x = element_blank()  # Remove x-axis title
  )
filepath <- glue("{man_image_path}/wholebrain_EIR_E4negative_4wayint_minage_age_sex_apoe4_MANUSCRIPT_NEWCOLORS.pdf")
ggsave(filepath, plot = plot1, width = 6, height = 3.5625) 
plot1


# Generate the plot for E4+ without y-axis title and legend
plot2 <- plot_model(wholebrain_lme4_E4_binary, type = "pred", terms = c("years_since_first_scan", "sex", "E4_binary[\u03B54+]"), axis.labels = "") 
color_levels <- levels(plot2$data$group_col)
color_mapping <- setNames(c(female_color, male_color), color_levels)
plot2 <- plot2 + 
  labs(
    title = NULL,
    x = "Years Since First Scan",
    y = NULL
  ) +
  scale_color_manual(
    name = "Sex",  # Change this to your desired legend title
    values = color_mapping  # Use the dynamically generated color mapping
  ) +
  scale_fill_manual(
    values = color_mapping  # Use the same color mapping for the confidence intervals
  ) +
  manuscript_theme +
  theme(
    axis.title.x = element_blank()  # Remove x-axis title
  )
filepath <- glue("{man_image_path}/wholebrain_EIR_E4positive_4wayint_minage_age_sex_apoe4_MANUSCRIPT_NEWCOLORS.pdf")
ggsave(filepath, plot = plot2, width = 6, height = 3.5625) #Special case: 1.14 x 1.92 inches. whole plot is 1.28 x 3.83 inches.
plot2




#===============================#
# * Plot histogram ----
#===============================#

### Plot histogram for all whole brain EIB
scaling_factor =3.125
midpoint <- mean(df_longitudinal$EIR_whole_brain)
nu_bins <- 13
breaks <- seq(1, 1.65, length.out = nu_bins + 1)
# Use ggplot_build to get the number of bins
gg_b <- ggplot_build(
  ggplot(df_longitudinal, aes(x = EIR_whole_brain)) +
    geom_histogram(breaks = breaks)
)
# Calculate the bin midpoints for gradient mapping
bin_midpoints <- seq(min(df_longitudinal$EIR_whole_brain), max(df_longitudinal$EIR_whole_brain), length.out = nu_bins)
# Create a custom gradient palette
custom_palette <- scales::gradient_n_pal(c(hist_low_color, hist_mid_color, hist_high_color))(scales::rescale(bin_midpoints))
# Generate the histogram with custom gradient colors
plot <- ggplot(df_longitudinal, aes(x = EIR_whole_brain)) +
  geom_histogram(breaks = breaks, aes(fill = ..x..), color = "black", show.legend = FALSE) +
  scale_fill_gradientn(
    colors = c(hist_low_color, hist_mid_color, hist_high_color),
    values = scales::rescale(c(min(df_longitudinal$EIR_whole_brain), mean(df_longitudinal$EIR_whole_brain), max(df_longitudinal$EIR_whole_brain)))
  ) +
  labs(
    x = "Whole Brain EIR",
    y = "Frequency"
  ) +
  manuscript_theme +
  theme(
    axis.title.y = element_blank()
  )
filepath <- glue("{man_image_path}/wholebrain_EIR_histogram_MANUSCRIPT_GRADIENT.pdf")
ggsave(filepath, plot = plot, width = 1.38*scaling_factor, height = 0.6*scaling_factor) #0.6 x 1.38 inches
plot




#===============================#
# * Plot Sphaghetti Plot for 3-way interaction  ----
#===============================#
p <- ggplot(df_longitudinal_agebysubjectcentered, aes(x = years_since_first_scan, y = EIR_whole_brain, group = subject_ID)) +
  geom_line(aes(color = sex), alpha = 0.5) + # Individual lines
  geom_smooth(aes(group = sex, color = sex), method = "lm", se = FALSE, size=2) + # LOESS or lm trend lines
  facet_grid(~E4_binary) + # Create 6 panels based on sex and age_category
  theme_minimal() + # Use a minimal theme
  labs(x = "Years Since First Scan", y = "Whole Brain EIR", color = "sex") # Labels
print(p)







#-#######################-#
# Model Asessment: glmmTMB ----
#-#######################-#

model_assessment_dir <- glue("{results_dir}/model_assessment")
if (!file.exists(model_assessment_dir)) {dir.create(model_assessment_dir, recursive = TRUE)}

#Define model
mod <- glmmTMB(EIR_whole_brain ~ scale(age_at_first_scan) * scale(years_since_first_scan) * sex * E4_binary + (1 | subject_ID), data = df_longitudinal_agebysubjectcentered)



#---------------------------------------------------#
# * Residuals against possible PREDICTORS ----
#---------------------------------------------------#
residuals_vs_predictors_dir <- glue("{model_assessment_dir}/residuals_vs_predictors")
if (!file.exists(residuals_vs_predictors_dir)) {dir.create(residuals_vs_predictors_dir, recursive = TRUE)}

plot <- residual_plots_predictors(data = df_longitudinal_agebysubjectcentered, modelTMB = mod, predictors = c(var_num, var_fac), y_axis_font_size=8) 
ggsave(glue("{residuals_vs_predictors_dir}/residuals_vs_predictors_wholebrain_APOEnpscore_all.pdf"), plot = plot, width = 10, height = 6, dpi = 300)

# NOTE: For this check, we recommend that var_num and var_fac ALSO contain available variable(s) 
# that are currently NOT part of the model. 

# What should I look for? 
# >> Points should homogeneously scatter without obvious pattern.

# Resolving violations:
# >> Seek a distribution family and/or link function that better capture the observed residual distribution.
# >> The model may lack an important covariate (shown here: covariates that show residual patterns).
#    or an informative interaction term (-> 6.1.3).
# >> The model may lack relevant non-linear (polynomial) terms (shown here: curvilinear patterns in residuals across a given covariate).
# >> The model could be over- or underdispersed (-> 6.2.1).
# >> The model may suffer from zero-inflation (-> 6.2.2).
# >> Add a dispersion formula to the model to explicitly integrate heterogeneous variance across predictor values / levels.



#------------------------------------------------------------------------------#  
# * Residuals against possible TWO-WAY INTERACTIONS AMONG PREDICTORS ----
#------------------------------------------------------------------------------# 
residuals_vs_2waypredictors_dir <- glue("{model_assessment_dir}/residuals_vs_2WayPredictors")
if (!file.exists(residuals_vs_2waypredictors_dir)) {dir.create(residuals_vs_2waypredictors_dir, recursive = TRUE)}

# Graphical displays of all 2-way predictor-response relationships.
savepath <- glue("{residuals_vs_2waypredictors_dir}/residuals_vs_2WayPredictors_wholebrain_APOEnpscore_all.pdf")
residual_plots_interactions(data = df_longitudinal_agebysubjectcentered, modelTMB = mod, predictors = c(var_num, var_fac), savepath = savepath) 

# NOTE-1: In numeric vs. numeric display panels, smoothing lines are added as a visual aid to detect non-linear relationships.
# NOTE-2: Up to 9 plots are shown in the plot window. 
#         If > 9 plots are present, they are saved in your working directory as "Two-ways_interactions.pdf".

# What should I look for?   
# >> Points should homogeneously scatter without obvious pattern. 

# Resolving violations:
# >> Consider adding the relevant interaction term to the model if meaningful. 




#=================================#
# * Autocorrelation checks ---- 
#=================================#
post_pred_checks_dir <- glue("{results_dir}/posterior_predictive_checks")
if (!file.exists(post_pred_checks_dir)) {dir.create(post_pred_checks_dir, recursive = TRUE)}

autocorrelation_dir <- glue("{post_pred_checks_dir}/autocorrelation")
if (!file.exists(autocorrelation_dir)) {dir.create(autocorrelation_dir, recursive = TRUE)}

# We check for spatial and temporal correlation patterns in model residuals (= autocorrelation) 
# using standardised semivariograms.
autocor_check(data = df_longitudinal_agebysubjectcentered, 
              modelTMB = mod,
              variable = var_time,       # Checking for TIME? => use var_time  | Checking for SPACE? => use var_space.
              grouping = var_time_groups, # (optional) - Grouping variable for multiple TIME series, only.
              # maxlag = NA,             # (optional) - Sets the maximum distance between pairs of observations to be inspected for autocorrelation. 
              n.sim = 500)               # Number of simulations, set to >2000 for final model assessment.

# What should I look for?
# >> Ideally, observed standardised semivariances should be WITHIN the pattern of permuted standardised 
#    semivariances extracted from the model.
# >> Autocorrelation typically shows when observed standardised semivariance falls clearly below 1 at 
#    shorter distances (= towards the left).

# Treating temporal/spatial autocorrelation? 
# >> Add an autocorrelation structure to the model formulation (see article section 5.2.2).  
# >> More detail on possible autocorrelation structures: 
#    https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html.








#-#######################-#
# Get summary statistics ----
#-#######################-#

#===============================#
# * All subjects (df) ----
#===============================#

# age by e4
age_by_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(
    mean.value = mean(age, na.rm=T),
    sd.value = sd(age, na.rm=T)
  )

# sex by e4
sex_by_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary, sex) %>%
  summarize(
    count = n_distinct(subject_ID)
  )

# number of scans by e4
scans_by_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary, subject_ID) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, count) %>%
  summarize(count_of_scan_counts = n())

# Number of e4
num_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(count = n()) 

# Number of subjects by e4
num_subjs_by_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(num_subjects = n_distinct(subject_ID)) 

# Number of amyloid scans
num_amyloid_scans_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(amyloid_binary) %>%
  summarize(count = n()) 

# Number of amyloid subjects (>1 scan)
num_amyloid_subj_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  summarize(any_positive = any(amyloid_binary == "Amyloid +")) %>%
  summarize(num_amyloid_positive_subjects = sum(any_positive),
            num_amyloid_negative_subjects = sum(!any_positive))

# Number of CDR
num_CDR_scans_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(CDR) %>%
  summarize(count = n()) 

# CDR by e4
num_CDR_scans_by_e4_table <- df_longitudinal_agebysubjectcentered %>%
  group_by(E4_binary, CDR) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, CDR, count) 


### Save tables
sink(file.path(glue("{results_dir}/wholebrain_binary_summary_statistics_tables.txt")))

cat("## Age by e4 table\n\n")
print(age_by_e4_table)
cat("## sex by e4 table\n\n")
print(sex_by_e4_table)
cat("\n\n## Number of scans by e4 table\n\n")
print(scans_by_e4_table)
cat("\n\n## Number of e4 table\n\n")
print(num_e4_table)
cat("\n\n## Number of subjects by e4 table\n\n")
print(num_subjs_by_e4_table)
cat("\n\n## Number of amyloid scans table\n\n")
print(num_amyloid_scans_table)
cat("\n\n## Number of amyloid subjects table\n\n")
print(num_amyloid_subj_table)
cat("\n\n## Number of CDR scans table\n\n")
print(num_CDR_scans_table)
cat("\n\n## CDR by e4 table\n\n")
print(num_CDR_scans_by_e4_table)

sink()



#===============================#
# * PET longituidnal (df_longitudinal_PET) ----
#===============================#

# age by e4
age_by_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(
    mean.value = mean(age, na.rm=T),
    sd.value = sd(age, na.rm=T)
  )

# sex by e4
sex_by_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary, sex) %>%
  summarize(
    count = n_distinct(subject_ID)
  )

# number of scans by e4
scans_by_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary, subject_ID) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, count) %>%
  summarize(count_of_scan_counts = n())

# Number of e4
num_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(count = n()) 

# Number of subjects by e4
num_subjs_by_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(num_subjects = n_distinct(subject_ID)) 

# Number of amyloid scans
num_amyloid_scans_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(amyloid_binary) %>%
  summarize(count = n()) 

# Number of amyloid subjects (>1 scan)
num_amyloid_subj_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  summarize(any_positive = any(amyloid_binary == "Amyloid +")) %>%
  summarize(num_amyloid_positive_subjects = sum(any_positive),
            num_amyloid_negative_subjects = sum(!any_positive))

# Number of CDR
num_CDR_scans_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(CDR) %>%
  summarize(count = n()) 

# CDR by e4
num_CDR_scans_by_e4_table <- df_longitudinal_PET_agebysubjectcentered %>%
  group_by(E4_binary, CDR) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, CDR, count) 


### Save tables
sink(file.path(glue("{results_dir}/wholebrain_PETlongitudinal_summary_statistics_tables.txt")))

cat("## Age by e4 table\n\n")
print(age_by_e4_table)
cat("## sex by e4 table\n\n")
print(sex_by_e4_table)
cat("\n\n## Number of scans by e4 table\n\n")
print(scans_by_e4_table)
cat("\n\n## Number of e4 table\n\n")
print(num_e4_table)
cat("\n\n## Number of subjects by e4 table\n\n")
print(num_subjs_by_e4_table)
cat("\n\n## Number of amyloid scans table\n\n")
print(num_amyloid_scans_table)
cat("\n\n## Number of amyloid subjects table\n\n")
print(num_amyloid_subj_table)
cat("\n\n## Number of CDR scans table\n\n")
print(num_CDR_scans_table)
cat("\n\n## CDR by e4 table\n\n")
print(num_CDR_scans_by_e4_table)

sink()



#===============================#
# * PET longituidnal (df_longitudinal_PET) ----
#===============================#

# age by e4
age_by_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(
    mean.value = mean(age, na.rm=T),
    sd.value = sd(age, na.rm=T)
  )

# sex by e4
sex_by_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary, sex) %>%
  summarize(
    count = n_distinct(subject_ID)
  )

# number of scans by e4
scans_by_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary, subject_ID) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, count) %>%
  summarize(count_of_scan_counts = n())

# Number of e4
num_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(count = n()) 

# Number of subjects by e4
num_subjs_by_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary) %>%
  summarize(num_subjects = n_distinct(subject_ID)) 

# Number of amyloid scans
num_amyloid_scans_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(amyloid_binary) %>%
  summarize(count = n()) 

# Number of amyloid subjects (>1 scan)
num_amyloid_subj_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(subject_ID) %>%
  summarize(any_positive = any(amyloid_binary == "Amyloid +")) %>%
  summarize(num_amyloid_positive_subjects = sum(any_positive),
            num_amyloid_negative_subjects = sum(!any_positive))

# Number of CDR
num_CDR_scans_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(CDR) %>%
  summarize(count = n()) 

# CDR by e4
num_CDR_scans_by_e4_table <- df_longitudinal_amyloidnegative_agebysubjectcentered %>%
  group_by(E4_binary, CDR) %>%
  summarize(count = n()) %>%
  group_by(E4_binary, CDR, count) 


### Save tables
sink(file.path(glue("{results_dir}/wholebrain_amyloidnegative_summary_statistics_tables.txt")))

cat("## Age by e4 table\n\n")
print(age_by_e4_table)
cat("## sex by e4 table\n\n")
print(sex_by_e4_table)
cat("\n\n## Number of scans by e4 table\n\n")
print(scans_by_e4_table)
cat("\n\n## Number of e4 table\n\n")
print(num_e4_table)
cat("\n\n## Number of subjects by e4 table\n\n")
print(num_subjs_by_e4_table)
cat("\n\n## Number of amyloid scans table\n\n")
print(num_amyloid_scans_table)
cat("\n\n## Number of amyloid subjects table\n\n")
print(num_amyloid_subj_table)
cat("\n\n## Number of CDR scans table\n\n")
print(num_CDR_scans_table)
cat("\n\n## CDR by e4 table\n\n")
print(num_CDR_scans_by_e4_table)

sink()




