slopes <- function(data, exclude_time = 0, window_duration = 60) {
  # Ensure the time column is numeric
  data$time <- as.numeric(data$time)

  # Filter the data to exclude the initial time period
  data <- data[data$time >= exclude_time, ]

  # Calculate the slope of MO2_blank vs time
  blank_lm <- lm(MO2_blank ~ time, data = data)
  slope_blank <- coef(blank_lm)[2]

  # Function to calculate slope, r-squared, and starting time in a rolling window
  roll_lm <- function(y, x, window_size) {
    slopes <- c()
    r_squared <- c()
    start_times <- c()

    n <- length(y)

    if (n < window_size) {
      stop("Window size is larger than the data length.")
    }

    for (i in 1:(n - window_size + 1)) {
      y_window <- y[i:(i + window_size - 1)]
      x_window <- x[i:(i + window_size - 1)]

      if (length(y_window) != window_size || length(x_window) != window_size) {
        next  # Skip this iteration if window sizes are incorrect
      }

      lm_fit <- lm(y_window ~ x_window)
      slopes <- c(slopes, coef(lm_fit)[2])
      r_squared <- c(r_squared, summary(lm_fit)$r.squared)
      start_times <- c(start_times, x[i])  # Store the starting time of the window
    }

    return(data.frame(slopes, r_squared, start_times))
  }

  # Determine the number of data points in the rolling window
  time_interval <- mean(diff(data$time))  # Assuming approximately constant time interval
  window_size <- ceiling(window_duration / time_interval)

  # Initialize a list to store results
  results_list <- list()

  # Loop over each MO2_fish column and calculate rolling slopes, r-squared, and start times
  fish_columns <- c("MO2_fish1", "MO2_fish2", "MO2_fish3", "MO2_fish4")
  for (fish_col in fish_columns) {
    if (!fish_col %in% names(data)) {
      warning(paste(fish_col, "is not in the data"))
      next
    }

    roll_results <- roll_lm(data[[fish_col]], data$time, window_size)

    # Adjust slopes by subtracting the slope of MO2_blank
    adjusted_slopes <- roll_results$slopes - slope_blank

    # Find the minimum adjusted slope
    min_adjusted_slope <- min(adjusted_slopes)
    min_slope_idx <- which.min(adjusted_slopes)
    min_slope_r_squared <- roll_results$r_squared[min_slope_idx]
    min_slope_start_time <- roll_results$start_times[min_slope_idx]
    mean_r_squared <- mean(roll_results$r_squared, na.rm = TRUE)

    # Calculate the percentage of the fish slope that the blank slope represents
    percentage_of_fish <- as.numeric(slope_blank / min_adjusted_slope) * 100

    results_list[[fish_col]] <- list(
      "Adjusted Slopes" = adjusted_slopes,
      "Starting Times" = roll_results$start_times,
      "Minimum Adjusted Slope" = min_adjusted_slope,
      "R-Squared for Minimum Slope" = min_slope_r_squared,
      "Start Time for Minimum Slope" = min_slope_start_time,
      "Mean R-Squared for all windows" = mean_r_squared,
      "Percentage of Fish Slope (%)" = percentage_of_fish
    )
  }

  return(results_list)
}
