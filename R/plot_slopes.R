plot_slopes <- function(data, exclude_time = 0, window_duration = 60) {
  # Ensure the time column is numeric
  if (!"time" %in% names(data)) {
    stop("The data frame must contain a 'time' column.")
  }
  data$time <- as.numeric(data$time)

  # Filter the data to exclude the initial time period
  data <- data[data$time >= exclude_time, ]

  # Check if the required columns are present
  required_columns <- c("MO2_blank", "MO2_fish1", "MO2_fish2", "MO2_fish3", "MO2_fish4")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop("The data frame is missing the following required columns: ", paste(missing_columns, collapse = ", "))
  }

  # Calculate the slope of MO2_blank vs time
  blank_lm <- lm(MO2_blank ~ time, data = data)
  slope_blank <- coef(blank_lm)[2]

  # Function to calculate slope, r-squared, and starting time in a rolling window
  roll_lm <- function(y, x, window_size) {
    n <- length(y)
    slopes <- numeric(n - window_size + 1)
    r_squared <- numeric(n - window_size + 1)
    start_times <- numeric(n - window_size + 1)

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
      slopes[i] <- coef(lm_fit)[2]
      r_squared[i] <- summary(lm_fit)$r.squared
      start_times[i] <- x[i]  # Store the starting time of the window
    }

    return(data.frame(slopes, r_squared, start_times))
  }

  # Determine the number of data points in the rolling window
  time_interval <- mean(diff(data$time))  # Assuming approximately constant time interval
  window_size <- ceiling(window_duration / time_interval)

  # Initialize a list to store results
  results_list <- list()
  slope_plot_list <- list()
  r2_plot_list <- list()
  all_r_squared <- c()  # Collect all r-squared values to determine global limits

  # Loop over each MO2_fish column and calculate rolling slopes, r-squared, and start times
  fish_columns <- c("MO2_fish1", "MO2_fish2", "MO2_fish3", "MO2_fish4")
  for (i in seq_along(fish_columns)) {
    fish_col <- fish_columns[i]
    if (!fish_col %in% names(data)) {
      warning(paste(fish_col, "is not in the data"))
      next
    }

    roll_results <- roll_lm(data[[fish_col]], data$time, window_size)

    # Adjust slopes by subtracting the slope of MO2_blank
    adjusted_slopes <- roll_results$slopes - slope_blank

    slope_colors = c("#4E79A7", "#919C4C", "#FD8F24", "#C03728")

    # Collect all r-squared values for global scaling
    all_r_squared <- c(all_r_squared, roll_results$r_squared)

    # Plot slopes and r-squared values against start times
    plot_data <- data.frame(
      Start_Times = roll_results$start_times,
      Slopes = adjusted_slopes,
      R_Squared = roll_results$r_squared
    )

    p_slope <- ggplot(plot_data, aes(x = Start_Times, y = Slopes)) +
      geom_line(color = slope_colors[i]) +
      labs(title = paste(fish_col, "Slopes vs Time"),
           x = "Time (s)",
           y = "Adjusted Slope (mg O2/L/s") +
      theme_minimal()

    slope_plot_list[[fish_col]] <- p_slope
    r2_plot_list[[fish_col]] <- plot_data
  }

  # Determine global limits for r-squared values
  r_squared_limits <- range(all_r_squared, na.rm = TRUE)

  # Plot r-squared values with global color scaling
  for (fish_col in fish_columns) {
    if (!is.null(r2_plot_list[[fish_col]])) {
      plot_data <- r2_plot_list[[fish_col]]
      p_r2 <- ggplot(plot_data, aes(x = Start_Times, y = R_Squared)) +
        geom_line(aes(color = R_Squared)) +
        scale_color_gradient(low = "red", high = "green", limits = r_squared_limits) +
        labs(title = paste(fish_col, "R-Squared vs Time"),
             x = "Time (s)",
             y = "R-Squared") +
        theme_minimal()
      r2_plot_list[[fish_col]] <- p_r2
    }
  }

  # Arrange slope plots in a grid and display
  slope_plots <- do.call(grid.arrange, c(slope_plot_list, ncol = 2))
  print(slope_plots)

  # Arrange r2 plots in a grid and display
  r2_plots <- do.call(grid.arrange, c(r2_plot_list, ncol = 2))
  print(r2_plots)
}
