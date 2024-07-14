window_effect <- function(data, exclude_time = 0, window_durations = seq(30, 600, 30), r_squared_threshold = 0.95) {
  results <- list()
  fish_columns <- c("MO2_fish1", "MO2_fish2", "MO2_fish3", "MO2_fish4")
  min_window_durations <- list()
  all_r_squared <- c()  # To collect all r-squared values for consistent scaling

  for (fish_col in fish_columns) {
    min_window_durations[[fish_col]] <- NA  # Initialize with NA in case no window duration meets the threshold
  }

  for (window_duration in window_durations) {
    cat("Analyzing window duration:", window_duration, "seconds\n")
    slopes_results <- slopes(data, exclude_time, window_duration)

    for (fish_col in fish_columns) {
      if (is.null(slopes_results[[fish_col]])) next
      if (!fish_col %in% names(results)) {
        results[[fish_col]] <- data.frame(Window_Duration = integer(),
                                          Minimum_Slope = numeric(),
                                          Mean_R_Squared = numeric())
      }
      results[[fish_col]] <- rbind(results[[fish_col]],
                                   data.frame(Window_Duration = window_duration,
                                              Minimum_Slope = slopes_results[[fish_col]]$`Minimum Adjusted Slope`,
                                              Mean_R_Squared = slopes_results[[fish_col]]$`Mean R-Squared for all windows`))

      mean_r_squared <- slopes_results[[fish_col]]$`Mean R-Squared for all windows`
      all_r_squared <- c(all_r_squared, mean_r_squared)

      if (!is.na(mean_r_squared) && mean_r_squared >= r_squared_threshold && is.na(min_window_durations[[fish_col]])) {
        min_window_durations[[fish_col]] <- window_duration
      }
    }
  }

  # Define common r-squared limits
  r_squared_limits <- range(all_r_squared, na.rm = TRUE)

  # Plotting the results
  plots <- list()
  for (fish_col in fish_columns) {
    if (!fish_col %in% names(results)) next
    p <- ggplot(results[[fish_col]], aes(x = Window_Duration, y = Minimum_Slope)) +
      geom_point(aes(color = Mean_R_Squared)) +
      scale_color_gradient(low = "red", high = "green", limits = r_squared_limits) +
      labs(title = paste(fish_col, "Window Analysis"),
           x = "Window Duration (seconds)",
           y = "Minimum Adjusted Slope") +
      theme_minimal()
    plots[[fish_col]] <- p
  }

  # Arrange plots in a grid
  grid.arrange(grobs = plots, ncol = 2)

  # Print minimum window durations for the specified R-squared threshold
  cat("\nMinimum window durations for R-squared threshold of", r_squared_threshold, ":\n")
  print(min_window_durations)
}
