


#=========================================================#
# Non-plotting functions ----
#=========================================================#
#' Calculate the Root Mean Square of Successive Differences (RMSSD)
#'
#' @param measurements A numeric vector of measurements. The vector can contain
#'   missing values (NA), which will be removed before calculations.
#'
#' @return A single numeric value representing the RMSSD of the input measurements.
#'   If the input vector has less than two non-missing values, the function returns NA.
#'
calculate_RMSSD <- function(measurements) {
  if (length(measurements) > 1) {
    measurements <- na.omit(measurements) # Calculate successive differences
    diff_values <- diff(measurements) # Calculate the squares of these differences
    squared_diffs <- diff_values^2 # Calculate the mean of the squared differences
    mean_squared_diffs = mean(squared_diffs, na.rm=TRUE) # Return the square root of the mean squared differences

    sqrt(mean_squared_diffs)
  } else {
    NA  # Return NA if there's not enough data to calculate RMSSD
  }
}





#=========================================================#
# Plotting Functions and Settings ----
#=========================================================#

manuscript_theme <- theme(
  text = element_text(family = "Arial"),
  plot.title = element_text(size = 20, face = "bold"),
  plot.subtitle = element_text(size = 22),
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  strip.text = element_text(size = 20, family = "Arial", face = "bold"),
  #strip.background = element_rect(fill = "lightblue", color = "grey", size = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  axis.line = element_line(color = "black", size = 0.5) 
)



