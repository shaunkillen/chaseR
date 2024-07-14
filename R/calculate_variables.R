calculate_variables <- function(data, tubing_diameters, tubing_lengths, V_chambers, fish_masses, exclude_time = 0, window_duration = 60, smr_values) {
  # Calculate slopes and MMR values
  slopes_results <- slopes(data, exclude_time, window_duration)

  fish_columns <- c("MO2_fish1", "MO2_fish2", "MO2_fish3", "MO2_fish4")
  results <- list()
  plots_MMR <- list()
  plots_EPOC <- list()
  summary_list <- list()
  MMR_values_list <- list()
  epoc_results <- list()

  colors <- c("#4E79A7", "#919C4C", "#FD8F24", "#C03728")  # Colors for the plots

  for (i in seq_along(fish_columns)) {
    fish_col <- fish_columns[i]
    if (!is.null(slopes_results[[fish_col]])) {
      adjusted_slopes <- slopes_results[[fish_col]]$`Adjusted Slopes`
      start_times <- slopes_results[[fish_col]]$`Starting Times`
      tubing_diameter <- tubing_diameters[i]
      tubing_length <- tubing_lengths[i]
      V_chamber <- V_chambers[i]
      fish_mass <- fish_masses[i]

      # Calculate volume in liters
      tubing_volume_liters <- (3.14159265359 * ((tubing_diameter / 2) ^ 2) * tubing_length) / 1000
      V_chamber_liters <- V_chamber / 1000
      fish_volume_liters <- fish_mass / 1000

      total_volume_liters <- tubing_volume_liters + V_chamber_liters - fish_volume_liters
      total_volume_ml <- total_volume_liters * 1000  # Convert back to mL

      # Calculate MMR
      MMR_values <- (-1 * adjusted_slopes) * 3600 * total_volume_liters
      max_MMR <- max(MMR_values)
      max_MMR_time <- start_times[which.max(MMR_values)]
      max_MMR_r_squared <- slopes_results[[fish_col]]$r_squared[which.max(MMR_values)]
      cat("Maximum Metabolic Rate for", fish_col, ":", max_MMR, "mg O2/h\n")
      cat("Total Volume for", fish_col, ":", total_volume_ml, "mL\n")

      # Store MMR results
      results[[fish_col]] <- list(
        "Maximum Metabolic Rate" = max_MMR,
        "R-Squared for MMR" = max_MMR_r_squared,
        "Time for MMR" = max_MMR_time,
        "Total Volume" = total_volume_ml,
        "Units" = list(MMR = "mg O2/h", Volume = "mL"),
        "MMR Values" = MMR_values,
        "Start Times" = start_times
      )

      # Define the model function
      model_func <- function(par, t) {
        par[1] * exp(-par[2] * t) + par[3]
      }

      # Fit the non-linear model
      start_params <- c(a = max(MMR_values), b = 0.01, c = min(MMR_values))
      nls_fit <- nlsLM(MMR_values ~ a * exp(-b * Time) + c,
                       start = start_params,
                       data = data.frame(Time = start_times, MMR_values = MMR_values))
      coef_fit <- coef(nls_fit)

      # Store non-linear model results
      results[[fish_col]]$`Non-linear Model Coefficients` <- coef_fit
      results[[fish_col]]$`Non-linear Model` <- nls_fit

      # Create a plot for MMR values over time with the fitted model
      plot_data <- data.frame(Time = start_times, MMR = MMR_values)
      fit_data <- data.frame(Time = seq(min(start_times), max(start_times), length.out = 100))
      fit_data$MMR <- predict(nls_fit, newdata = fit_data)

      p_MMR <- ggplot(plot_data, aes(x = Time, y = MMR)) +
        geom_point(size=0.5, color = colors[i]) +
        geom_line(data = fit_data, aes(x = Time, y = MMR), color = colors[i], linetype = "dashed") +
        geom_hline(yintercept = max_MMR, linetype = "dashed", color = colors[i]) +  # Add dashed horizontal line
        annotate("text", x = max_MMR_time, y = max_MMR, label = paste("MMR =", round(max_MMR, 2), "mg O2/h"),
                 vjust = -0.6, hjust = -0.5, color = colors[i]) +
        expand_limits(y = max(MMR_values) * 1.1) +  # Expand y-axis to ensure label fits
        labs(title = paste("MO2 vs Time for", fish_col),
             x = "Time (s)",
             y = "MO2 (mg O2/h)") +
        theme_minimal()
      plots_MMR[[fish_col]] <- p_MMR

      # Store summary results
      summary_list[[fish_col]] <- list(
        "Maximum Metabolic Rate (mg O2/h)" = max_MMR,
        "R-Squared for Maximum MMR" = max_MMR_r_squared,
        "Time for MMR (s)" = max_MMR_time,
        "Exponential Rate of Recovery (1/s)" = coef_fit["b"],
        "Asymptote (mg O2/h)" = coef_fit["c"]
      )

      # Calculate EPOC if SMR values are provided
      if (!is.null(smr_values)) {
        smr <- smr_values[i]

        # Adjust MMR values by subtracting the SMR
        adjusted_MMR_values <- MMR_values - smr

        # Calculate the area under the curve using the trapezoid rule
        auc <- sum(diff(start_times) * (head(adjusted_MMR_values, -1) + tail(adjusted_MMR_values, -1)) / 2)

        # Store the EPOC result
        epoc_results[[fish_col]] <- auc

        # Add EPOC plot
        plot_data_EPOC <- data.frame(Time = start_times, MMR = MMR_values)
        p_EPOC <- ggplot(plot_data_EPOC, aes(x = Time, y = MMR)) +
          geom_line(color = colors[i]) +
          geom_hline(yintercept = smr, linetype = "dashed", color = "black") +
          labs(title = paste("EPOC for", fish_col),
               x = "Time (s)",
               y = "MO2 (mg O2/h)") +
          theme_classic()

        # Add trapezoids to the plot
        for (j in 1:(length(start_times) - 1)) {
          trapezoid_data <- data.frame(
            x = c(start_times[j], start_times[j + 1], start_times[j + 1], start_times[j]),
            y = c(MMR_values[j], MMR_values[j + 1], smr, smr)
          )
          p_EPOC <- p_EPOC + geom_polygon(data = trapezoid_data, aes(x = x, y = y), fill = colors[i], alpha = 0.2)
        }

        plots_EPOC[[fish_col]] <- p_EPOC

        # Add EPOC to summary results
        summary_list[[fish_col]]$EPOC <- auc
      }
    }
  }

  # Arrange MMR plots in a grid
  grid.arrange(grobs = plots_MMR, ncol = 2)

  # Arrange EPOC plots in a grid if available
  if (length(plots_EPOC) > 0) {
    grid.arrange(grobs = plots_EPOC, ncol = 2)
  }

  # Print the summary results
  for (fish_col in fish_columns) {
    if (!is.null(summary_list[[fish_col]])) {
      cat("\nSummary for", fish_col, ":\n")
      cat("Maximum Metabolic Rate (mg O2/h):", summary_list[[fish_col]]$`Maximum Metabolic Rate (mg O2/h)`, "\n")
      cat("R-Squared for MMR:", summary_list[[fish_col]]$`R-Squared for MMR`, "\n")
      cat("Time for MMR (s):", summary_list[[fish_col]]$`Time for MMR (s)`, "\n")
      cat("Exponential Rate of Recovery (1/s):", summary_list[[fish_col]]$`Exponential Rate of Recovery (1/s)`, "\n")
      cat("Asymptote (mg O2/h):", summary_list[[fish_col]]$`Asymptote (mg O2/h)`, "\n")
      if (!is.null(summary_list[[fish_col]]$EPOC)) {
        cat("EPOC (mg O2):", summary_list[[fish_col]]$EPOC, "\n")
      }
    }
  }
}
