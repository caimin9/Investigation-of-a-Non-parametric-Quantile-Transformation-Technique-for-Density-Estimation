# Load required packages
library(dplyr)
library(MASS)
library(parallel)      # For parallel processing
library(doParallel)    # For foreach parallel backend
library(foreach)       # For foreach loops
options(pillar.sigfig = 4)

# Set up parallel backend
num_cores = detectCores() - 1  # Leave one core free for system processes
cl = makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary functions and libraries to the workers
clusterEvalQ(cl, {
  library(dplyr)
  library(MASS)
})

run_combined_simulation = function(n, seed, shape, scale, initial_df) {
  set.seed(seed)
  
  # Generate data using input shape and scale parameters
  Y = rweibull(n, shape = shape, scale = scale)
  
  
  upper_bound = qweibull(0.9999, shape = shape, scale = scale)
  # Add a buffer to ensure we capture the full distribution
  upper_bound = upper_bound 
  
  # Create adaptive grid
  y_grid = seq(0, upper_bound, length.out = 1000)
  true_density = dweibull(y_grid, shape = shape, scale = scale)
  
  #------------------------------------------------------------------------------
  # 1. Common setup
  #------------------------------------------------------------------------------
  probs = ppoints(n)
  normal_q = qnorm(probs)
  
  # Helper functions
  spline_weights = function(quantiles, n, density_values) {
    return(quantiles * (1 - quantiles) / (n * (density_values^2)))
  }
  
  # Trapezoidal rule for numerical integration
  trapz = function(x, y) {
    # Approximate the integral of y w.r.t. x using the trapezoidal rule.
    sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(x)]) / 2)
  }
  
  # Transformation function
  transformation = function(spline_obj, x_vals, f_x_in) { 
    dy_dx = predict(spline_obj, x = x_vals, deriv = 1)$y
    f_y = f_x_in / abs(dy_dx)
    return(f_y)
  }
  
  #------------------------------------------------------------------------------
  # 2. Iterative Method
  #------------------------------------------------------------------------------
  sample_q = quantile(Y, probs)
  f_x = dnorm(normal_q)
  
  
  # 4) Iterative procedure
  iterative_qde = function(x_vals, y_vals, x_dist, probs, n, df, true_density, num_iter = 30) {
    # Start with uniform weights = 1
    weights = spline_weights(probs, n, x_dist)
    weights_matrix = matrix(NA, nrow = n, ncol = num_iter+1)
    error_vec = numeric(num_iter)
    weights_matrix[,1] = weights
    for (iter in seq_len(num_iter)) {
      spline_fit = smooth.spline(x = x_vals, y = y_vals, w = 1 / weights, df = df)
      f_y_new = transformation(spline_fit, x_vals, f_x)
      updated_weights = spline_weights(probs, n, f_y_new) # Update weights
      weights = updated_weights
      weights_matrix[, iter+1] = weights
      
      # Calculate ISE at this iteration (using the global y_grid)
      y_pred_iter = spline_fit$y
      f_y_interp_iter = approx(y_pred_iter, f_y_new, xout = y_grid, rule = 2)$y
      err_squared = (f_y_interp_iter - true_density)^2
      error_vec[iter] = trapz(y_grid, err_squared)
    }
    
    # Return the final spline and final f_y
    list(spline_fit = spline_fit, final_f_y  = f_y_new, weights = weights, weights_matrix = weights_matrix,
         error_per_iter = error_vec)
  } # end of function
  
  # Run the iterative QDE with specified initial df
  res_iter = iterative_qde(x_vals = normal_q, y_vals = sample_q, x_dist = f_x, probs = probs, 
                           n = n, true_density = true_density, df = initial_df, num_iter = 30)
  
  ss = res_iter$spline_fit
  iterative_f_y = res_iter$final_f_y
  final_weights = res_iter$weights_matrix[, which.min(res_iter$error_per_iter)]
  
  
  
  # ISE function for optimisation
  calculate_ISE = function(df_value) {
    spline_fit = smooth.spline(normal_q, sample_q, w = 1 / final_weights, df = df_value)
    y_pred     = spline_fit$y
    f_y_pred   = transformation(spline_fit, normal_q, f_x)
    f_y_interp = approx(y_pred, f_y_pred, xout = y_grid, rule = 2)$y
    
    squared_diff = (f_y_interp - true_density)^2
    
    # Integrate the squared difference over y via the trapezoidal rule
    ise_val = trapz(y_grid, squared_diff)
    return(ise_val)
  }
  optim_df_iter = optim(par = 4, fn = calculate_ISE, method = "Brent", lower = 1, upper = 30)
  
  optim_df = optim(par = 4, fn = calculate_ISE, method = "Brent", lower = 1, upper = 30)
  optimal_df = optim_df$par
  
  # Final density estimate (need it for plots)
  spline_final = smooth.spline(normal_q, sample_q, w = 1 / final_weights, df = optimal_df)
  y_pred_final = spline_final$y
  f_y_trans = transformation(spline_final, normal_q, f_x)
  fY_iter_final = approx(y_pred_final, f_y_trans, xout = y_grid, rule = 2)$y
  
  ise_direct = optim_df$value
  
  #------------------------------------------------------------------------------
  # 3. Direct Method (Non-iterative)
  #------------------------------------------------------------------------------
  sample_q = quantile(Y, probs = probs)
  f_x = dnorm(normal_q)
  weights_direct = spline_weights(probs, n, f_x)
  
  # ISE function for optimisation
  calculate_ISE_direct = function(df_value) {
    spline_fit = smooth.spline(normal_q, sample_q, w = 1 / weights_direct, df = df_value)
    y_pred     = spline_fit$y
    f_y_pred   = transformation(spline_fit, normal_q, f_x)
    f_y_interp = approx(y_pred, f_y_pred, xout = y_grid, rule = 2)$y
    
    squared_diff = (f_y_interp - true_density)^2
    
    # Integrate the squared difference over y via the trapezoidal rule
    ise_val = trapz(y_grid, squared_diff)
    return(ise_val)
  }
  # Optimise df for direct method
  optim_df_direct = optim(par = 4, fn = calculate_ISE_direct, method = "Brent", lower = 1, upper = 30)
  
  # Final direct method density
  spline_direct = smooth.spline(normal_q, sample_q, w = 1/weights_direct, df = optim_df_direct$par)
  y_pred_direct = predict(spline_direct, x = normal_q)$y
  fY_est_direct = transformation(spline_direct, normal_q, f_x)
  fY_direct_final = approx(y_pred_direct, fY_est_direct, xout = y_grid, rule = 2)$y
  
  #------------------------------------------------------------------------------
  # 4. KDE Method
  #------------------------------------------------------------------------------
  calculate_ISE_kde = function(bw_value) {
    kde_fit = density(Y, kernel = 'gaussian', bw = bw_value,n =n)
    kde_pdf = approx(kde_fit$x, kde_fit$y, xout = y_grid, rule = 2)$y
    squared_diff = (kde_pdf - true_density)^2
    return(trapz(y_grid, squared_diff))
  }
  
  optim_bw = optim(par = 1.0, fn = calculate_ISE_kde, method = "Brent", lower = 0.01, upper = 6.0)
  
  # Final KDE
  kde_final = density(Y, bw = optim_bw$par, kernel = "gaussian",n=n)
  kde_pdf = approx(kde_final$x, kde_final$y, xout = y_grid, rule = 2)$y
  
  
  #------------------------------------------------------------------------------
  # 4. FITDISTR Method
  #------------------------------------------------------------------------------
  fit_mle = fitdistr(Y, densfun = "weibull")
  mle_shape = fit_mle$estimate["shape"]
  mle_rate  = fit_mle$estimate["scale"]
  
  # Construct the fitted MLE weibull PDF on the y_grid
  fd_pdf = dweibull(y_grid, shape = mle_shape, scale = mle_rate)
  fd_pdf[is.infinite(fd_pdf)] = 1e-30
  
  
  #------------------------------------------------------------------------------
  # Calculate final MSEs
  #------------------------------------------------------------------------------
  iter_sq_error = (fY_iter_final - true_density)^2
  iter_ise = trapz(y_grid,iter_sq_error)
  direct_sq_error = (fY_direct_final - true_density)^2
  direct_ise = trapz(y_grid,direct_sq_error)
  fd_sq_error = (fd_pdf - true_density)^2
  fd_ise = trapz(y_grid,fd_sq_error)
  kde_squared_error = (kde_pdf - true_density)^2
  kde_ise = trapz(y_grid,kde_squared_error)
  
  
  return(list(
    n = n,
    seed = seed,
    shape = shape,
    scale = scale,
    initial_df = initial_df,
    iter_ise = iter_ise,
    direct_ise = direct_ise,
    fd_ise = fd_ise,
    kde_ise = kde_ise,
    optimal_df_iter = optim_df_iter$par,
    optimal_df_direct = optim_df_direct$par,
    optimal_bw = optim_bw$par,
    densities = list(
      y_grid = y_grid,
      true_density = true_density,
      iter_density = fY_iter_final,
      direct_density = fY_direct_final,
      fd_density = fd_pdf,
      kde_density = kde_pdf
    )
  ))
}


# Run simulations in parallel
n_values = 100 * 2^(0:4)
seeds = 1:100
weibull_params = list(
  list(shape = 1.2, scale = 2)
  #list(shape = 1.5, scale = 3),
  #list(shape = 2.0, scale = 1.5),
  #list(shape = 2.5, scale = 1.2),
  #list(shape = 3, scale = 1.8)   
)

# Define different df values to test
initial_df_values = c(2,4,6,8,10,12,14,16,18,20)

# Create expanded results storage
results_all = list()

# Modified simulation loop using parallel processing
for(gp in weibull_params) {
  shape_val = gp$shape
  scale_val = gp$scale
  
  for(initial_df in initial_df_values) {
    cat(sprintf("\nRunning simulations for shape=%g, scale=%g, initial_df=%g\n", 
                shape_val, scale_val, initial_df))
    
    results_combined = data.frame()
    density_examples = list()
    
    # Parallelize across n and seed combinations
    # Create all combinations of parameters
    param_grid = expand.grid(n = n_values, seed = seeds)
    
    # Run parallel foreach loop
    parallel_results = foreach(i = 1:nrow(param_grid), .combine = rbind, 
                                .packages = c("dplyr", "MASS")) %dopar% {
                                  n = param_grid$n[i]
                                  seed = param_grid$seed[i]
                                  
                                  # For progress tracking (will print from each worker)
                                  cat(sprintf("Running simulation: n=%d, seed=%d on worker %d\n", 
                                              n, seed, Sys.getpid()))
                                  
                                  # Run the simulation
                                  sim_result = run_combined_simulation(n, seed, shape_val, scale_val, initial_df)
                                  
                                  # Return results as data frame
                                  data.frame(
                                    shape = shape_val,
                                    scale = scale_val,
                                    initial_df = initial_df,
                                    n = sim_result$n,
                                    seed = sim_result$seed,
                                    iter_ise = sim_result$iter_ise,
                                    direct_ise = sim_result$direct_ise,
                                    fd_ise = sim_result$fd_ise,
                                    kde_ise = sim_result$kde_ise,
                                    df_iter = sim_result$optimal_df_iter,
                                    df_direct = sim_result$optimal_df_direct,
                                    bw = sim_result$optimal_bw,
                                    # Add a flag for seed==1 to identify which results should save densities
                                    save_density = (seed == 1)
                                  )
                                }
    
    # Combine results
    results_combined = parallel_results
    
    # Extract density examples separately (only for seed==1)
    for(n in n_values) {
      # Run one more time for seed 1 to get densities
      # (alternatively, you could modify the parallel code to return densities,
      # but this keeps the memory usage lower during parallel execution)
      sim_result = run_combined_simulation(n, 1, shape_val, scale_val, initial_df)
      density_examples[[as.character(n)]] = sim_result$densities
    }
    
    # Store results with df in key
    key = sprintf("shape_%g_scale_%g_df_%g", shape_val, scale_val, initial_df)
    results_all[[key]] = list(
      results = results_combined,
      densities = density_examples
    )
  }
}


stopCluster(cl)

#########################################################################################################
################ RESULTS ANALYSIS
#########################################################################################################


#results_all = readRDS("Weibull_results_all.rds") #Do this if data already saved


# Save results
# 1. Create boxplots for each sample size
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
for(n in n_values) {
  subset_data = results_combined[results_combined$n == n,]
  boxplot(subset_data$iter_ise, subset_data$direct_ise, subset_data$fd_ise,subset_data$kde_ise,
          names = c("Iterative", "Direct","FD", "KDE"),
          main = sprintf("n = %d", n),
          ylab = "MSE",
          col = c("lightblue", "lightgreen", "pink"))
}

# 2. Create density comparison plots
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
for(n in n_values) {
  dens = density_examples[[as.character(n)]]
  plot(dens$y_grid, dens$true_density, type = "l", col = "black", lwd = 2,
       main = sprintf("Density Estimates (n = %d)", n),
       xlab = "y", ylab = "Density")
  lines(dens$y_grid, dens$iter_density, col = "blue", lwd = 2, lty = 2)
  lines(dens$y_grid, dens$direct_density, col = "red", lwd = 2, lty = 3)
  lines(dens$y_grid, dens$kde_density, col = "green", lwd = 2, lty = 4)
  lines(dens$y_grid, dens$fd_density, col = "yellow", lwd = 2, lty = 5)
  
  legend("topright", 
         legend = c("True", "Iterative", "Direct", "KDE", "Fit Dist"),
         col = c("black", "blue", "red", "green", 'yellow'),
         lty = 1:4, lwd = 2)
}
################################################
## Change dot types for better plotting in B&W

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

# Define point characters for each estimator
pch_vals = c(NA, 16, 17, 15, 18)  # solid circle, triangle, square, diamond

for(n in n_values) {
  dens = density_examples[[as.character(n)]]
  
  # Plot the true density as a plain line
  plot(dens$y_grid, dens$true_density, type = "l", col = "black", lwd = 2,
       main = sprintf("Density Estimates (n = %d)", n),
       xlab = "y", ylab = "Density")
  
  # Add each estimator with line and points (same color, matching symbols)
  lines(dens$y_grid, dens$iter_density, col = "blue", lwd = 2)
  points(dens$y_grid, dens$iter_density, col = "blue", pch = 16)
  
  lines(dens$y_grid, dens$kde_density, col = "green", lwd = 2)
  points(dens$y_grid, dens$kde_density, col = "green", pch = 17)
  
  lines(dens$y_grid, dens$direct_density, col = "red", lwd = 2)
  points(dens$y_grid, dens$direct_density, col = "red", pch = 15)
  
  lines(dens$y_grid, dens$fd_density, col = "orange", lwd = 2)
  points(dens$y_grid, dens$fd_density, col = "orange", pch = 18)
  
  # Add legend (matching symbols and line colors)
  legend("topright", 
         legend = c("True", "Iterative", "KDE", "Normal Weights", "Fit Dist"),
         col = c("black", "blue", "green", "red", "orange"),
         lty = c(1, 1, 1, 1, 1),
         pch = c(NA, 16, 17, 15, 18),
         lwd = 2)
}











########################################

# 3. Summary statistics
summary_stats = results_combined %>%
  group_by(n) %>%
  summarise(
    iter_mean = mean(iter_ise),iter_sd = sd(iter_ise),direct_mean = mean(direct_ise), 
    direct_sd = sd(direct_ise),kde_mean = mean(kde_ise), kde_sd = sd(kde_ise), 
    fd_mean = mean(fd_ise), fd_sd = sd(fd_ise))

print("Summary Statistics:")
print(summary_stats)

# 4. Calculate relative performance
for(n in n_values) {
  subset_data = results_combined[results_combined$n == n,]
  cat(sprintf("\nResults for n = %d:\n", n))
  cat("Iterative vs KDE wins:", mean(subset_data$iter_ise < subset_data$kde_ise), "\n")
  cat("Direct vs KDE wins:", mean(subset_data$direct_ise < subset_data$kde_ise), "\n")
  cat("Iterative vs Direct wins:", mean(subset_data$iter_ise <= subset_data$direct_ise), "\n")
  
  # Mean percentage improvements
  iter_improvement = mean((subset_data$kde_ise - subset_data$iter_ise) / subset_data$kde_ise * 100)
  direct_improvement = mean((subset_data$kde_ise - subset_data$direct_ise) / subset_data$kde_ise * 100)
  
  cat(sprintf("Mean improvement Iterative over KDE: %.2f%%\n", iter_improvement))
  cat(sprintf("Mean improvement Direct over KDE: %.2f%%\n", direct_improvement))
}

# 5. Plot optimal parameters
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
boxplot(df_iter ~ n, data = results_combined,
        main = "Optimal df (Iterative)",
        xlab = "Sample Size", ylab = "Degrees of Freedom")
boxplot(df_direct ~ n, data = results_combined,
        main = "Optimal df (Direct)",
        xlab = "Sample Size", ylab = "Degrees of Freedom")
boxplot(bw ~ n, data = results_combined,
        main = "Optimal Bandwidth (KDE)",
        xlab = "Sample Size", ylab = "Bandwidth")




options(scipen =999)
  
# Analysis across all parameter combinations
for(gp in weibull_params) {
  shape_val = gp$shape
  scale_val = gp$scale
  cat(sprintf("\nAnalysis for weibull(%g, %g):\n", shape_val, scale_val))
  
  # Collect results for all df values for this shape/scale combination
  df_results = data.frame()
  for(initial_df in initial_df_values) {
    key = sprintf("shape_%g_scale_%g_df_%g", shape_val, scale_val, initial_df)
    df_results = rbind(df_results, results_all[[key]]$results)
  }
  
  # Analyze performance by initial df value
  df_summary = df_results %>%
    group_by(initial_df, n) %>%
    summarise(
      mean_ise = mean(iter_ise),
      sd_ise = sd(iter_ise)
    )
  
  print("Summary by initial df value:")
  print(df_summary)
  
  # Create visualisation
  par(mfrow = c(1, 1))



  # Define a better color palette (more distinguishable than rainbow)
  colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
  
  # Set up the plot with improved parameters
  par(mar = c(5, 5, 4, 2))  # Increase margins for labels
  plot(NULL, xlim = range(initial_df_values), 
       ylim = range(df_summary$mean_ise),
       main = bquote("Mean ISE vs Initial df Weibull("*.(shape_val)*","*.(scale_val)*")"),
       xlab = "Initial df", 
       ylab = "",
       cex.lab = 1.2,      # Larger axis labels
       cex.axis = 1.1,     # Larger axis numbers
       cex.main = 1.3,     # Larger title
       las = 1)            # Horizontal y-axis labels
  
  # Add grid for better readability
  grid(nx = NULL, ny = NULL, lty = "dotted", col = "gray80")
  
  # Plot lines with improved aesthetics
  for(n_val in n_values) {
    subset_data = df_summary[df_summary$n == n_val,]
    lines(subset_data$initial_df, subset_data$mean_ise, 
          col = colors[which(n_values == n_val)],
          type = "o",                  # Points and lines
          pch = 16,                    # Solid circle points
          lwd = 2,                     # Thicker lines
          cex = 1.2)                   # Larger points
  }
  
  # Improved legend
  legend("topright", 
         legend = paste("n =", format(n_values, big.mark=",")), # Add thousands separator
         col = colors,
         lty = 1,
         pch = 16,                     # Match points style
         pt.cex = 1.2,                 # Match point size
         lwd = 2,                      # Match line width
         bty = 'n',                    # No box
         inset = c(0.02, 0.1),
         bg = "white",                 # White background
         box.col = "white",            # White box edge
         cex = 1.1)                    # Larger text
  
  title(ylab = "Mean ISE", line = 4, cex.lab = 1.2)  # Increase line value to move label down
  
  # Add box around plot
  box()
}

#####################################
# Define a better color palette and point styles
colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
point_styles = c(16, 17, 15, 18, 3, 8)  # Solid circle, triangle, square, diamond, plus, star

# Set up the plot
par(mar = c(5, 5, 4, 8))
plot(NULL, xlim = range(initial_df_values), 
     ylim = range(log(df_summary$mean_ise)),
     main = bquote("Mean ISE vs Initial df: Weibull("*.(shape_val)*","*.(scale_val)*")"),
     xlab = "Initial df", 
     ylab = "",
     cex.lab = 1.2,
     cex.axis = 1.1,
     cex.main = 1.3,
     las = 1)
grid(nx = NULL, ny = NULL, lty = "dotted", col = "gray80")

# Loop through n-values with different point characters
for(i in seq_along(n_values)) {
  n_val = n_values[i]
  subset_data = df_summary[df_summary$n == n_val, ]
  lines(subset_data$initial_df, log(subset_data$mean_ise), 
        col = colors[i],
        type = "o",
        pch = point_styles[i],
        lwd = 2,
        cex = 1.2)
}



title(ylab = "Log(MISE)", line = 4, cex.lab = 1.2)
box()

legend("topright", 
       inset = c(-0.3, 0),  # Push legend outside
       xpd = TRUE,           # Allow drawing outside plot region
       legend = paste("n =", format(n_values, big.mark=",")),
       col = colors,
       pch = point_styles,
       lty = 1,
       lwd = 2,
       pt.cex = 1.2,
       cex = 1.1,
       bty = 'n')
#####################################



results_df_choice = results_combined %>%
  group_by(initial_df) %>%
  summarise(
    mean_iter_ise = mean(iter_ise),
    sd_iter_ise   = sd(iter_ise),
    n_cases       = n()
  )
results_df_choice


# Assuming each element of results_all has a component "results"
all_results = do.call(rbind, lapply(results_all, function(x) x$results))

# Now group by initial_df to summarise iter_ise across all simulations
results_df_choice = all_results %>%
  group_by(initial_df) %>%
  summarise(
    mean_iter_ise = median(iter_ise),
    sd_iter_ise   = sd(iter_ise),
    n_cases       = n()
  )

results_df_choice






################################################################################
############ Line Plots of ISE vs DF
################################################################################
par(mfrow = c(2, 2))

for(gp in weibull_params) {
  shape_val = gp$shape
  scale_val = gp$scale
  cat(sprintf("\nAnalysis for weibull(%g, %g):\n", shape_val, scale_val))
  
  # Collect results for all df values for this shape/scale combination
  df_results = data.frame()
  for(initial_df in initial_df_values) {
    key = sprintf("shape_%g_scale_%g_df_%g", shape_val, scale_val, initial_df)
    df_results = rbind(df_results, results_all[[key]]$results)
  }
  
  # Analyse performance by initial df value
  df_summary = df_results %>%
    group_by(initial_df, n) %>%
    summarise(
      mean_ise = mean(iter_ise),
      sd_ise = sd(iter_ise)
    )
  
  print("Summary by initial df value:")
  print(df_summary)
  
  # Create visualization
  
  
  # Line plot of mean ISE vs initial df for different n
  plot(NULL, xlim = range(initial_df_values), 
       ylim = range(df_summary$mean_ise),
       main = sprintf("Mean ISE vs Initial df\nGamma(%g, %g)", shape_val, scale_val),
       xlab = "Initial df", ylab = "Mean ISE")
  
  for(n_val in n_values) {
    subset_data = df_summary[df_summary$n == n_val,]
    lines(subset_data$initial_df, subset_data$mean_ise, 
          col = rainbow(length(n_values))[which(n_values == n_val)],
          type = "b")
  }
  legend("topright", legend = paste("n =", n_values),
         col = rainbow(length(n_values)), lty = 1, bty ='n')
}





#################################################################################
##### PLOT HOW ISE CHANGES VS INITIAL CHOICE OF DF
#################################################################################
# Combine all results into one data frame
df_all = do.call(rbind, lapply(names(results_all), function(nm) {
  df = results_all[[nm]]$results
  # Create a distribution label using the first row's shape and scale values
  df$dist_label = paste0("weibull(", df$shape[1], ",", df$scale[1], ")")
  return(df)
}))

# Remove rows where initial_df is 100
df_all = subset(df_all, initial_df != 100)

# Aggregate the data: compute mean iterative ISE for each distribution, n, and initial_df
agg_data = aggregate(iter_ise ~ dist_label + n + initial_df,
                     data = df_all, FUN = mean)

# Get the unique distribution labels (e.g. "weibull(6,6)", "weibull(1,4)", etc.)
dist_levels = unique(agg_data$dist_label)

# Set up the layout for the plots (here 1 row and 2 columns)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))  # adjust margins if needed

# Loop over each distribution and plot mean iterative ISE vs. initial_df on a log scale
for(dl in dist_levels) {
  sub = agg_data[agg_data$dist_label == dl, ]
  # Get unique sample sizes within this distribution
  sample_sizes = unique(sub$n)
  
  # Set up an empty plot with a log-scaled y-axis.
  plot(NA, NA, xlim = range(sub$initial_df), ylim = range(sub$iter_ise),
       log = "y", xlab = "Initial Degrees of Freedom", 
       ylab = "Mean Iterative ISE", main = dl, type = "n")
  
  # Choose colours for different sample sizes
  cols = rainbow(length(sample_sizes))
  
  # Loop over each sample size and add a line (with points)
  for(i in seq_along(sample_sizes)) {
    ss = sample_sizes[i]
    dat = sub[sub$n == ss, ]
    dat = dat[order(dat$initial_df), ]  # ensure points are in order
    lines(dat$initial_df, dat$iter_ise, type = "b", col = cols[i], pch = 16)
  }
  
  # Compute positions for the legend: place it near the right edge and vertically centred.
  xlims = range(sub$initial_df)
  ylims = range(sub$iter_ise)
  # x position: a bit to the left of the maximum x value
  x_pos = xlims[2] - 0.2 * diff(xlims)
  # y position: geometric mean of the y-range (since log scale)
  y_pos = exp((log(ylims[1]) + log(ylims[2]))/2.05)
  
  legend(x = x_pos, y = y_pos, legend = sample_sizes, col = cols, lty = 1, pch = 16,
         title = "Sample Size (n)", cex = 1, bty = "n")
}



par(mfrow = c(1,1))
# Loop over all parameter combinations in shape_params, scale_params, and initial_df_values
for(gp in weibull_params) {
  shape_val = gp$shape
  scale_val = gp$scale
  cat(sprintf("\nAnalysis for Weibull(%g, %g):\n", shape_val, scale_val))
  
  # Collect results for all df values for this shape/scale combination
  df_results = data.frame()
  for(initial_df in initial_df_values) {
    key = sprintf("shape_%g_scale_%g_df_%g", shape_val, scale_val, initial_df)
    df_results = rbind(df_results, results_all[[key]]$results)
  }
  
  # Analyze performance by initial df value
  df_summary = df_results %>%
    group_by(initial_df, n) %>%
    summarise(
      mean_ise = mean(iter_ise),
      sd_ise = sd(iter_ise)
    )
  
  print("Summary by initial df value:")
  print(df_summary)
  
  # Define a better color palette
  colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
  
  # Set up the plot with improved parameters
  par(mar = c(5, 6, 4, 8))  # Increased right margin for legend
  plot(NULL, xlim = range(initial_df_values), 
       ylim = range(log(df_summary$mean_ise)),
       main = bquote("Mean ISE vs Initial df Weibull("*.(shape_val)*","*.(scale_val)*")"),
       xlab = "Initial df", 
       ylab = "",
       cex.lab = 1.2,      # Larger axis labels
       cex.axis = 1.1,     # Larger axis numbers
       cex.main = 1.3,     # Larger title
       las = 1,            # Horizontal y-axis labels
       mgp = c(3, 1, 0),   # Adjust margin line for axis labels and ticks
       log = "x")          # Log scale for x-axis
  
  # Add grid with appropriate spacing
  abline(h = seq(min(log(df_summary$mean_ise)), 
                 max(log(df_summary$mean_ise)), 
                 length.out = 10), 
         v = initial_df_values,
         col = "gray80", lty = "dotted")
  
  # Plot lines with improved aesthetics
  for(n_val in n_values) {
    subset_data = df_summary[df_summary$n == n_val,]
    lines(subset_data$initial_df, log(subset_data$mean_ise), 
          col = colors[which(n_values == n_val)],
          type = "o",                  # Points and lines
          pch = 16,                    # Solid circle points
          lwd = 2,                     # Thicker lines
          cex = 1.2)                   # Larger points
  }
  
  # Improved legend outside plot
  legend("topright", 
         legend = paste("n =", format(n_values, big.mark=",")), # Add thousands separator
         col = colors,
         lty = 1,
         pch = 16,                     # Match points style
         pt.cex = 1.2,                 # Match point size
         lwd = 2,                      # Match line width
         bty = 'n',                    # No box
         inset = c(-0.325, 0),           # Move legend outside
         xpd = TRUE,                   # Allow plotting outside
         cex = 1.1)                    # Larger text
  
  # Add y-axis label with adjusted position
  title(ylab = "log(Mean ISE)", line = 4.5, cex.lab = 1.2)
  
  # Add box around plot
  box()
}








plot(summary_stats$n,log(summary_stats$iter_mean))

#################################################################################### 
############################ Statistics bY RESULTS 
####################################################################################
all_results_df = do.call(rbind, lapply(results_all, function(x) x$results))

# 2. Summarise 
## Do 4 Groupings (order matters on how you write them)
# by n
# by shape, scale , n
# by nothing
# by shape,scale
overall_summary = all_results_df %>%
  group_by(shape,scale,n) %>%
  summarise(
    mean_iter_ise   = mean(iter_ise),
    mean_direct_ise = mean(direct_ise),
    mean_kde_ise    = mean(kde_ise),
    mean_fd_ise     = mean(fd_ise),
    sd_iter_ise   = sd(iter_ise),
    sd_direct_ise = sd(direct_ise),
    sd_kde_ise    = sd(kde_ise),
    sd_fd_ise     = sd(fd_ise),
    
    .groups         = "drop"   # So dplyr create nested groups
  )

# 3. Print or examine the summary
print(n=1000, overall_summary)

names(all_results_df)



################################################################################
######## Relative Efficiency
##############################################################################

summary_table = all_results_df %>%
  group_by(shape, scale,n) %>%
  summarize(
    mean_iter_vs_direct = mean(iter_ise/direct_ise),
    mean_direct_vs_kde = mean(direct_ise/kde_ise),
    mean_iter_vs_kde = mean(iter_ise/kde_ise),
  ) %>%
  ungroup()
print(summary_table, n =100)



####### PLOTS###################


# First, create the relative efficiency columns in your original dataset
all_results_df = all_results_df %>%
  mutate(
    iter_vs_direct = iter_ise/direct_ise,
    direct_vs_kde = direct_ise/kde_ise,
    iter_vs_kde = iter_ise/kde_ise
  )

# For each distribution (shape/scale combination)
for(shp in unique(all_results_df$shape)) {
  for(scl in unique(all_results_df$scale[all_results_df$shape == shp])) {
    
    # Filter data for this specific distribution
    dist_data = all_results_df %>% 
      filter(shape == shp, scale == scl)
    
    # Calculate mean relative efficiency for each sample size
    mean_data = dist_data %>%
      group_by(n) %>%
      summarize(
        mean_iter_vs_direct = mean(iter_vs_direct),
        mean_direct_vs_kde = mean(direct_vs_kde),
        mean_iter_vs_kde = mean(iter_vs_kde)
      )
    
    # Create a line plot with all three metrics
    # Set up the plot with explicit x-axis values
    plot(mean_data$n, mean_data$mean_iter_vs_direct,
         type = "n", 
         log = "x",  # Log scale for x-axis
         xlim = range(mean_data$n),
         ylim = range(c(mean_data$mean_iter_vs_direct, 
                        mean_data$mean_direct_vs_kde, 
                        mean_data$mean_iter_vs_kde)),
         xlab = "Sample Size", 
         ylab = "Mean Relative Efficiency",
         main = sprintf("Relative Efficiency for Weibull(%g, %g)", shp, scl),
         xaxt = "n") # Suppress default x-axis
    
    # Create custom x-axis with actual values
    axis(1, at = unique(mean_data$n), labels = unique(mean_data$n))
    
    # Add reference line at y = 1
    abline(h = 1, lty = 2, col = "gray")
    
    # Add the three lines
    lines(mean_data$n, mean_data$mean_iter_vs_direct, 
          type = "o", col = "blue", pch = 16)
    lines(mean_data$n, mean_data$mean_direct_vs_kde, 
          type = "o", col = "green", pch = 17)
    lines(mean_data$n, mean_data$mean_iter_vs_kde, 
          type = "o", col = "red", pch = 18)
    
    # Add legend
    legend("topright", 
           legend = c("Iterative vs Direct", "Direct vs KDE", "Iterative vs KDE"),
           col = c("blue", "green", "red"), 
           pch = c(16, 17, 18), 
           lty = 1,
           bty = "n")
  }
}



#######################################################################
########## Iter vs KDE
###########################################################################
comparison_data = all_results_df %>%
  group_by(shape, scale, n) %>%
  summarize(
    mean_iter_vs_kde = mean(iter_vs_kde),
    .groups = "drop"
  )

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_iter_vs_kde),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = "Iterative Weights vs KDE Relative Efficiency",
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray")

# Define colors for different distributions
dist_colors = rainbow(length(unique(paste0(comparison_data$shape, "-", comparison_data$scale))))
dist_counter = 1

# Add lines for each distribution
for(shp in unique(comparison_data$shape)) {
  for(scl in unique(comparison_data$scale[comparison_data$shape == shp])) {
    dist_data = comparison_data %>% 
      filter(shape == shp, scale == scl) %>%
      arrange(n) # Make sure data is ordered by sample size
    
    lines(dist_data$n, dist_data$mean_iter_vs_kde, 
          type = "o", 
          col = dist_colors[dist_counter], 
          pch = 15 + dist_counter)
    
    dist_counter = dist_counter + 1
  }
}

# Add legend
legend_labels = paste0("Weibull(", 
                        unique(comparison_data$shape), 
                        ",", 
                        unique(comparison_data$scale), 
                        ")")
legend("topright", 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       bty = "n")





#######################################################################
########## NW vs KDE
###########################################################################

comparison_data = all_results_df %>%
  group_by(shape, scale, n) %>%
  summarize(
    mean_direct_vs_kde = mean(direct_vs_kde),
    .groups = "drop"
  )

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_direct_vs_kde),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = "Normal Weights vs KDE Relative Efficiency",
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray")

# Define colors for different distributions
dist_colors = rainbow(length(unique(paste0(comparison_data$shape, "-", comparison_data$scale))))
dist_counter = 1

# Add lines for each distribution
for(shp in unique(comparison_data$shape)) {
  for(scl in unique(comparison_data$scale[comparison_data$shape == shp])) {
    dist_data = comparison_data %>% 
      filter(shape == shp, scale == scl) %>%
      arrange(n) # Make sure data is ordered by sample size
    
    lines(dist_data$n, dist_data$mean_direct_vs_kde, 
          type = "o", 
          col = dist_colors[dist_counter], 
          pch = 15 + dist_counter)
    
    dist_counter = dist_counter + 1
  }
}

# Add legend
legend_labels = paste0("Weibull(", 
                        unique(comparison_data$shape), 
                        ",", 
                        unique(comparison_data$scale), 
                        ")")
legend("topright", 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       bty = "n")





#######################################################################
########## For Iter vs direct
###########################################################################

comparison_data = all_results_df %>%
  group_by(shape, scale, n) %>%
  summarize(
    mean_iter_vs_direct = mean(iter_vs_direct),
    .groups = "drop"
  )

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_iter_vs_direct),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = "Iterative Weights vs Normal Weights Relative Efficiency",
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray")

# Define colors for different distributions
dist_colors = rainbow(length(unique(paste0(comparison_data$shape, "-", comparison_data$scale))))
dist_counter = 1

# Add lines for each distribution
for(shp in unique(comparison_data$shape)) {
  for(scl in unique(comparison_data$scale[comparison_data$shape == shp])) {
    dist_data = comparison_data %>% 
      filter(shape == shp, scale == scl) %>%
      arrange(n) # Make sure data is ordered by sample size
    
    lines(dist_data$n, dist_data$mean_iter_vs_direct, 
          type = "o", 
          col = dist_colors[dist_counter], 
          pch = 15 + dist_counter)
    
    dist_counter = dist_counter + 1
  }
}

# Add legend
legend_labels = paste0("Weibull(", 
                        unique(comparison_data$shape), 
                        ",", 
                        unique(comparison_data$scale), 
                        ")")
legend("topright", 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       bty = "n")










################################################################################
######## LOG(ISE) VS LOG(N)
##############################################################################


##############################################################################################################
##############                        ESTIMATION CHARCTERISTICS
##############################################################################################################


############################################################################################################
############### Log of iter ISE VS n with Larger Boxplots
############################################################################################################

par(mfrow = c(1,1))
weibull_params = unique(all_results_df[, c("shape", "scale")])
for (i in 1:nrow(weibull_params)) {
  
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Calculate mean for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(iter_ise)),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Calculate segment slopes between adjacent points
  segment_slopes = numeric(length(x_values) - 1)
  for (j in 1:(length(x_values) - 1)) {
    segment_slopes[j] = (y_values[j+1] - y_values[j]) / (x_values[j+1] - x_values[j])
  }
  
  # Same colour scheme as before
  box_color = "#B3D9FF"  # Light blue for boxplots
  line_color = "#0066CC"  # Darker blue for the line
  point_color = "#003366"  # Even darker blue for points
  ref_color = "#CC3366"    # Pink/red for reference line
  slope_color = "#228B22"  # Green for slope values
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$iter_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = substitute(paste("Log of Iterative ISE vs Sample Size - ", Weibull(s, t)), 
                         list(s = current_shape, t = current_scale)),
       xlab = "Sample Size (Log Scale) ", 
       ylab = "ISE (Log Scale)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  # Add grid before plotting data
  grid(lty = 1, col = "gray90")
  
  # Draw the axes with n_values instead of log(n_values)
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Add boxplots at each sample size with improved appearance
  for (n_val in n_values) {
    n_subset = subset(df_subset, n == n_val)
    boxplot(log(n_subset$iter_ise), add = TRUE, at = log(n_val),
            boxwex = 0.3, col = box_color, outline = TRUE,
            outcol = "gray40", outpch = 20, outcex = 0.5, 
            axes = FALSE, border = "#5588BB")
  }
  
  # Add the main line after boxplots
  lines(log(summary_stats$n), summary_stats$mean_log_ise,
        col = line_color, lwd = 2)
  
  # Add points with a different colour
  points(log(summary_stats$n), summary_stats$mean_log_ise,
         pch = 16, col = point_color, cex = 1.1)
  
  # Add a reference line showing n^(-4/5) convergence rate
  n_ref = n_values[1]
  ise_ref = exp(summary_stats$mean_log_ise[summary_stats$n == n_ref])
  theoretical_line = log(ise_ref) - 4/5 * (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = ref_color, lwd = 2, lty = 2)
  
  # Add segment slope annotations
  for (j in 1:(length(x_values) - 1)) {
    mid_x = (x_values[j] + x_values[j+1]) / 2
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 0.5  # Position above the line
    text(mid_x, mid_y, sprintf("%.2f", segment_slopes[j]), 
         cex = 0.8, col = slope_color)
  }
  
  # Change the text box position to top right
  text_box = par("usr")[2] - 0.3 * diff(par("usr")[1:2])
  text_box_y = par("usr")[4] - 0.25 * diff(par("usr")[3:4])
  
  # Create a better looking text box for the slope info
  rect(text_box - 1.5, 
       text_box_y - 0.25, 
       text_box + 1.5, 
       text_box_y + 0.25, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray70")
  
  # Add the slope value text
  text(text_box, text_box_y, 
       sprintf("Overall slope: %.3f (Theoretical: -0.800)", overall_slope),
       cex = 0.9, col = "black")
  
  # Add legend with no box (bty='n')
  legend("topright", 
         legend = c(bquote(Weibull(.(current_shape), .(current_scale))), 
                    expression(n^(-4/5)~"rate"),
                    "Segment slopes"),
         col = c(line_color, ref_color, slope_color),
         lwd = c(2, 2, NA),
         lty = c(1, 2, NA),
         pch = c(16, NA, NA),
         pt.cex = c(1.1, NA, NA),
         text.col = c("black", "black", slope_color),
         bty = 'n',
         cex = 1.0)
}



############################################################################################################
############### Log of kde ISE VS n with Larger Boxplots
############################################################################################################

par(mfrow = c(1,1))
weibull_params = unique(all_results_df[, c("shape", "scale")])
for (i in 1:nrow(weibull_params)) {
  
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Calculate mean for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(kde_ise)),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Calculate segment slopes between adjacent points
  segment_slopes = numeric(length(x_values) - 1)
  for (j in 1:(length(x_values) - 1)) {
    segment_slopes[j] = (y_values[j+1] - y_values[j]) / (x_values[j+1] - x_values[j])
  }
  
  # Same colour scheme as before
  box_color = "#B3D9FF"  # Light blue for boxplots
  line_color = "#0066CC"  # Darker blue for the line
  point_color = "#003366"  # Even darker blue for points
  ref_color = "#CC3366"    # Pink/red for reference line
  slope_color = "#228B22"  # Green for slope values
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$kde_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding-.5, max(y_range) + y_padding),
       main = substitute(paste("Log of KDE ISE vs Sample Size - ", Weibull(s, t)), 
                         list(s = current_shape, t = current_scale)),
       xlab = "Sample Size (Log Scale) ", 
       ylab = "ISE (Log Scale)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  
  # Add grid before plotting data
  grid(lty = 1, col = "gray90")
  
  # Draw the axes with n_values instead of log(n_values)
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Add boxplots at each sample size with improved appearance
  for (n_val in n_values) {
    n_subset = subset(df_subset, n == n_val)
    boxplot(log(n_subset$kde_ise), add = TRUE, at = log(n_val),
            boxwex = 0.3, col = box_color, outline = TRUE,
            outcol = "gray40", outpch = 20, outcex = 0.5, 
            axes = FALSE, border = "#5588BB")
  }
  
  # Add the main line after boxplots
  lines(log(summary_stats$n), summary_stats$mean_log_ise,
        col = line_color, lwd = 2)
  
  # Add points with a different colour
  points(log(summary_stats$n), summary_stats$mean_log_ise,
         pch = 16, col = point_color, cex = 1.1)
  
  # Add a reference line showing n^(-4/5) convergence rate
  n_ref = n_values[1]
  ise_ref = exp(summary_stats$mean_log_ise[summary_stats$n == n_ref])
  theoretical_line = log(ise_ref) - 4/5 * (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = ref_color, lwd = 2, lty = 2)
  
  # Add segment slope annotations
  for (j in 1:(length(x_values) - 1)) {
    mid_x = (x_values[j] + x_values[j+1]) / 2
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 0.5  # Position above the line
    text(mid_x, mid_y, sprintf("%.2f", segment_slopes[j]), 
         cex = 0.8, col = slope_color)
  }
  
  # Change the text box position to top right
  text_box = par("usr")[2] - 0.3 * diff(par("usr")[1:2])
  text_box_y = par("usr")[4] - 0.25 * diff(par("usr")[3:4])
  
  # Create a better looking text box for the slope info
  rect(text_box - 1.5, 
       text_box_y - 0.25, 
       text_box + 1.5, 
       text_box_y + 0.25, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray70")
  
  # Add the slope value text
  text(text_box, text_box_y, 
       sprintf("Overall slope: %.3f (Theoretical: -0.800)", overall_slope),
       cex = 0.9, col = "black")
  
  # Add legend with no box (bty='n')
  legend("topright", 
         legend = c(bquote(Weibull(.(current_shape), .(current_scale))), 
                    expression(n^(-4/5)~"rate"),
                    "Segment slopes"),
         col = c(line_color, ref_color, slope_color),
         lwd = c(2, 2, NA),
         lty = c(1, 2, NA),
         pch = c(16, NA, NA),
         pt.cex = c(1.1, NA, NA),
         text.col = c("black", "black", slope_color),
         bty = 'n',
         cex = 1.0)
}



############################################################################################################
############### Log of Direct ISE VS n with Larger Boxplots
############################################################################################################

par(mfrow = c(1,1))
weibull_params = unique(all_results_df[, c("shape", "scale")])
for (i in 1:nrow(weibull_params)) {
  
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Calculate mean for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(direct_ise)),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Calculate segment slopes between adjacent points
  segment_slopes = numeric(length(x_values) - 1)
  for (j in 1:(length(x_values) - 1)) {
    segment_slopes[j] = (y_values[j+1] - y_values[j]) / (x_values[j+1] - x_values[j])
  }
  
  
  
  
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$direct_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = substitute(paste("Log of NW ISE vs Sample Size - ", Weibull(s, t)), 
                         list(s = current_shape, t = current_scale)),
       xlab = "Sample Size (Log Scale) ", 
       ylab = "ISE (Log Scale)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  # Add grid before plotting data
  grid(lty = 1, col = "gray90")
  
  # Draw the axes with n_values instead of log(n_values)
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Add boxplots at each sample size with improved appearance
  for (n_val in n_values) {
    n_subset = subset(df_subset, n == n_val)
    boxplot(log(n_subset$direct_ise), add = TRUE, at = log(n_val),
            boxwex = 0.3, col = box_color, outline = TRUE,
            outcol = "gray40", outpch = 20, outcex = 0.5, 
            axes = FALSE, border = "#5588BB")  # Darker orange border
  }
  
  # Add the main line after boxplots
  lines(log(summary_stats$n), summary_stats$mean_log_ise,
        col = line_color, lwd = 2)
  
  # Add points with a different colour
  points(log(summary_stats$n), summary_stats$mean_log_ise,
         pch = 16, col = point_color, cex = 1.1)
  
  # Add a reference line showing n^(-4/5) convergence rate
  n_ref = n_values[1]
  ise_ref = exp(summary_stats$mean_log_ise[summary_stats$n == n_ref])
  theoretical_line = log(ise_ref) - 4/5 * (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = ref_color, lwd = 2, lty = 2)
  
  # Add segment slope annotations with increased vertical offset
  for (j in 1:(length(x_values) - 1)) {
    mid_x = (x_values[j] + x_values[j+1]) / 2
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 0.8  # Increased vertical offset
    text(mid_x, mid_y, sprintf("%.2f", segment_slopes[j]), 
         cex = 0.8, col = slope_color)
  }
  
  # Add text box for slope info
  text_box = par("usr")[2] - 0.3 * diff(par("usr")[1:2])
  text_box_y = par("usr")[4] - 0.25 * diff(par("usr")[3:4])
  
  # Create a better looking text box for the slope info
  rect(text_box - 1.5, 
       text_box_y - 0.25, 
       text_box + 1.5, 
       text_box_y + 0.25, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray70")
  
  # Add the slope value text with confidence intervals
  text(text_box, text_box_y, 
       sprintf("Overall slope: %.3f (Theoretical: -0.800)", overall_slope),
       cex = 0.9, col = "black")
  
  # Add legend with no box (bty='n')
  legend("topright", 
         legend = c(bquote(Weibull(.(current_shape), .(current_scale))), 
                    expression(n^(-4/5)~"rate"),
                    "Segment slopes"),
         col = c(line_color, ref_color, slope_color),
         lwd = c(2, 2, NA),
         lty = c(1, 2, NA),
         pch = c(16, NA, NA),
         pt.cex = c(1.1, NA, NA),
         text.col = c("black", "black", slope_color),
         bty = 'n',
         cex = 1.0)
}

############################################################################################################
############### Log of FD ISE VS n with Larger Boxplots
############################################################################################################

par(mfrow = c(1,1))
weibull_params = unique(all_results_df[, c("shape", "scale")])
for (i in 1:nrow(weibull_params)) {
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Calculate mean for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(fd_ise)),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Calculate segment slopes between adjacent points
  segment_slopes = numeric(length(x_values) - 1)
  for (j in 1:(length(x_values) - 1)) {
    segment_slopes[j] = (y_values[j+1] - y_values[j]) / (x_values[j+1] - x_values[j])
  }
  
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$fd_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = substitute(paste("Log of FD ISE vs Sample Size - ", Weibull(s, t)), 
                         list(s = current_shape, t = current_scale)),
       xlab = "Sample Size (Log Scale) ", 
       ylab = "ISE (Log Scale)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  # Add grid before plotting data
  grid(lty = 1, col = "gray90")
  
  # Draw the axes with n_values instead of log(n_values)
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Add boxplots at each sample size with improved appearance
  for (n_val in n_values) {
    n_subset = subset(df_subset, n == n_val)
    boxplot(log(n_subset$fd_ise), add = TRUE, at = log(n_val),
            boxwex = 0.3, col = box_color, outline = TRUE,
            outcol = "gray40", outpch = 20, outcex = 0.5, 
            axes = FALSE, border = "#5588BB")  # Darker purple border
  }
  
  # Add the main line after boxplots
  lines(log(summary_stats$n), summary_stats$mean_log_ise,
        col = line_color, lwd = 2)
  
  # Add points with a different colour
  points(log(summary_stats$n), summary_stats$mean_log_ise,
         pch = 16, col = point_color, cex = 1.1)
  
  # Add a reference line showing n^(-4/5) convergence rate
  n_ref = n_values[1]
  ise_ref = exp(summary_stats$mean_log_ise[summary_stats$n == n_ref])
  theoretical_line = log(ise_ref) - 4/5 * (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = ref_color, lwd = 2, lty = 2)
  
  # Add segment slope annotations with increased vertical offset
  for (j in 1:(length(x_values) - 1)) {
    mid_x = (x_values[j] + x_values[j+1]) / 2
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 1.6  # Increased vertical offset
    text(mid_x, mid_y, sprintf("%.2f", segment_slopes[j]), 
         cex = 0.8, col = slope_color)
  }
  
  # Add text box for slope info
  text_box = par("usr")[2] - 0.3 * diff(par("usr")[1:2])
  text_box_y = par("usr")[4] - 0.25 * diff(par("usr")[3:4])
  
  # Create a better looking text box for the slope info
  rect(text_box - 1.5, 
       text_box_y - 0.25, 
       text_box + 1.5, 
       text_box_y + 0.25, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray70")
  
  # Add the slope value text
  text(text_box, text_box_y, 
       sprintf("Overall slope: %.3f (Theoretical: -0.800)", overall_slope),
       cex = 0.9, col = "black")
  
  # Add legend with no box (bty='n')
  legend("topright", 
         legend = c(bquote(Weibull(.(current_shape), .(current_scale))), 
                    expression(n^(-4/5)~"rate"),
                    "Segment slopes"),
         col = c(line_color, ref_color, slope_color),
         lwd = c(2, 2, NA),
         lty = c(1, 2, NA),
         pch = c(16, NA, NA),
         pt.cex = c(1.1, NA, NA),
         text.col = c("black", "black", slope_color),
         bty = 'n',
         cex = 1.0)
}






























##############################################################################################################
################## ALL THE ESTIMATOR BOXPLOTS AT ONCE
##############################################################################################################
# Function to analyze convergence rates
analyze_convergence = function(df_subset, n_values, estimator_name = "IW", ise_column = "iter_ise") {
  # Calculate mean and other statistics for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(!!sym(ise_column))),
      median_log_ise = median(log(!!sym(ise_column))),
      sd_log_ise = sd(log(!!sym(ise_column))),
      q25_log_ise = quantile(log(!!sym(ise_column)), 0.25),
      q75_log_ise = quantile(log(!!sym(ise_column)), 0.75),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  conf_int = confint(lm_fit, "x_values", level = 0.95)
  
  # Calculate segment slopes between adjacent points
  segment_slopes = numeric(length(x_values) - 1)
  for (j in 1:(length(x_values) - 1)) {
    segment_slopes[j] = (y_values[j+1] - y_values[j]) / (x_values[j+1] - x_values[j])
  }
  
  # Calculate theoretical convergence rate for this estimator
  theoretical_rate = switch(estimator_name,
                             "IW" = -4/5,
                             "KDE" = -4/5,
                             "NW" = -4/5,
                             "FD" = -1,
                             -4/5)  # Default to -4/5
  
  # Return all analysis results
  return(list(
    summary_stats = summary_stats,
    x_values = x_values,
    y_values = y_values,
    overall_slope = overall_slope,
    conf_int = conf_int,
    segment_slopes = segment_slopes,
    theoretical_rate = theoretical_rate
  ))
}

# Function to create enhanced plots with cleaner visualization
create_enhanced_plot = function(df_subset, analysis_results, shape, scale, 
                                 n_values, estimator_name = "IW", ise_column = "iter_ise") {
  # Extract analysis results
  summary_stats = analysis_results$summary_stats
  x_values = analysis_results$x_values
  y_values = analysis_results$y_values
  overall_slope = analysis_results$overall_slope
  conf_int = analysis_results$conf_int
  segment_slopes = analysis_results$segment_slopes
  theoretical_rate = analysis_results$theoretical_rate
  
  # Get full estimator name
  full_estimator_name = switch(estimator_name,
                                "IW" = "Iterative Weights",
                                "NW" = "Normal Weights",
                                "KDE" = "KDE",
                                "FD" = "FD",
                                estimator_name)
  
  # Updated color scheme with more distinguishable colors
  box_color = switch(estimator_name,
                      "IW" = "#A1D8E6",    # Light blue
                      "NW" = "#FFB366",    # Light orange
                      "KDE" = "#B3E6B3",   # Light green
                      "FD" = "#E6B3CC",    # Light purple
                      "#A1C9F4")           # Default light blue
  
  line_color = switch(estimator_name,
                       "IW" = "#0066CC",    # Dark blue
                       "NW" = "#FF6600",    # Orange
                       "KDE" = "#009933",   # Dark green
                       "FD" = "#CC3366",    # Dark purple
                       "#0072B2")           # Default blue
  
  point_color = line_color
  ref_color = "#D55E00"   # Orange-red for reference line
  slope_color = "#228B22"  # Green for slope values
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset[[ise_column]]), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(x_values, y_values,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = paste0(full_estimator_name, ": Log ISE vs Sample Size - Weibull(", shape, ",", scale, ")"),
       xlab = "Sample Size (log scale)", 
       ylab = "log(ISE)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  # Add grid before plotting data
  grid(lty = 1, col = "gray90")
  
  # Draw the axes with n_values instead of log(n_values)
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Add boxplots at each sample size with improved appearance
  for (n_val in n_values) {
    n_subset = subset(df_subset, n == n_val)
    boxplot(log(n_subset[[ise_column]]), add = TRUE, at = log(n_val),
            boxwex = 0.3, col = box_color, outline = TRUE,
            outcol = "gray40", outpch = 20, outcex = 0.5, 
            axes = FALSE, border = line_color)
  }
  
  # Add the main line after boxplots
  lines(x_values, y_values,
        col = line_color, lwd = 2)
  
  # Add points with a different colour
  points(x_values, y_values,
         pch = 16, col = point_color, cex = 1.1)
  
  # Add a reference line showing theoretical convergence rate
  n_ref = n_values[1]
  ise_ref = exp(summary_stats$mean_log_ise[summary_stats$n == n_ref])
  theoretical_line = log(ise_ref) + theoretical_rate * (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = ref_color, lwd = 2, lty = 2)
  
  # Add segment slope annotations with increased vertical offset
  for (j in 1:(length(x_values) - 1)) {
    mid_x = (x_values[j] + x_values[j+1]) / 2
    # Increased vertical offset from 0.3 to 0.8
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 0.8  
    text(mid_x, mid_y, sprintf("%.2f", segment_slopes[j]), 
         cex = 0.8, col = slope_color)
  }
  
  # Add text box for slope info
  text_box_x = mean(par("usr")[1:2])
  text_box_y = par("usr")[3] + 0.15 * diff(par("usr")[3:4])
  
  # Create a cleaner text box
  rect(par("usr")[1] + 0.05 * diff(par("usr")[1:2]), 
       text_box_y - 0.4, 
       par("usr")[2] - 0.05 * diff(par("usr")[1:2]), 
       text_box_y + 0.4, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray80")
  
  # Add the slope value text with statistical information
  text(text_box_x, text_box_y, 
       sprintf("Empirical slope: %.3f (95%% CI: [%.3f, %.3f])", 
               overall_slope, conf_int[1], conf_int[2]),
       cex = 0.9, col = "black")
  text(text_box_x, text_box_y - 0.3, 
       sprintf("Theoretical rate: %.3f", theoretical_rate),
       cex = 0.9, col = "black")
  
  # Add legend
  legend("topright", 
         legend = c(paste0(estimator_name), 
                    paste0("Theoretical n^(", theoretical_rate, ")"),
                    "Segment slopes"),
         col = c(line_color, ref_color, slope_color),
         lwd = c(2, 2, NA),
         lty = c(1, 2, NA),
         pch = c(16, NA, NA),
         pt.cex = c(1.1, NA, NA),
         text.col = c("black", "black", slope_color),
         bty = 'n',
         cex = 1.0)
}

# Main analysis loop for each Weibull parameter combination
weibull_params = unique(all_results_df[, c("shape", "scale")])
for (i in 1:nrow(weibull_params)) {
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  # Get subset of results for current parameters
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Create a list to store results for all estimators
  all_estimator_results = list()
  
  # Set up plot grid
  par(mfrow = c(2, 2))
  
  # Analyze each estimator
  estimators = c("IW", "KDE", "NW", "FD")
  ise_columns = c("iter_ise", "kde_ise", "direct_ise", "fd_ise")
  
  for (e in 1:length(estimators)) {
    estimator = estimators[e]
    ise_col = ise_columns[e]
    
    # Analyze convergence
    analysis_results = analyze_convergence(df_subset, n_values, estimator, ise_col)
    all_estimator_results[[estimator]] = analysis_results
    
    # Create plot
    create_enhanced_plot(df_subset, analysis_results, current_shape, current_scale, 
                         n_values, estimator, ise_col)
  }
  
  # Create comparison plot
  par(mfrow = c(1, 1))
  
  # Get y-range for all estimators
  all_y_values = unlist(lapply(all_estimator_results, function(x) x$y_values))
  y_range = range(all_y_values)
  y_padding = 0.1 * diff(y_range)
  
  # Set up comparison plot
  plot(log(n_values), log(n_values), type = "n",
       xlim = range(log(n_values)),
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = paste0("Comparison of Estimators - Weibull(", current_shape, ",", current_scale, ")"),
       xlab = "Sample Size (log scale)", 
       ylab = "log(ISE)",
       cex.main = 1.4, 
       cex.lab = 1.2, 
       cex.axis = 1.1,
       axes = FALSE)
  
  # Add grid and axes
  grid(lty = 1, col = "gray90")
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  box()
  
  # Updated colors for better distinction
  colors = c("#0066CC",    # Dark blue (IW)
              
              "#009933",    # Dark green (NW)
              "#FF6600",# Orange (KDE)
              "#CC3366")    # Dark purple (FD)
  ltys = c(1, 2, 3, 4)
  
  for (e in 1:length(estimators)) {
    estimator = estimators[e]
    results = all_estimator_results[[estimator]]
    
    # Plot line and points
    lines(results$x_values, results$y_values, 
          col = colors[e], lwd = 2, lty = ltys[e])
    points(results$x_values, results$y_values,
           pch = 15 + e, col = colors[e], cex = 1.1)
  }
  
  # Add legend with full estimator names and slopes
  legend_text = sapply(1:length(estimators), function(e) {
    estimator = estimators[e]
    results = all_estimator_results[[estimator]]
    
    full_name = switch(estimator,
                        "IW" = "Iterative Weights",
                        "NW" = "Normal Weights",
                        "KDE" = "KDE",
                        "FD" = "FD",
                        estimator)
    
    sprintf("%s (slope: %.3f)", full_name, results$overall_slope)
  })
  
  legend("topright", legend = legend_text,
         col = colors, lwd = 2, lty = ltys, pch = 15 + 1:length(estimators),
         pt.cex = 1.1, bty = 'n', cex = 1.0)
}



####################################################################
######
##################################################################
# Function to generate a LaTeX table of results for the Weibull distribution
generate_latex_table_weibull = function(all_estimator_results, current_shape, current_scale) {
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{Convergence Rate Analysis for Weibull(", current_shape, ",", current_scale, ")}\n")
  cat("\\begin{tabular}{lccc}\n")
  cat("\\hline\n")
  cat("Estimator & Empirical Rate & Theoretical Rate & 95\\% CI \\\\\n")
  cat("\\hline\n")
  
  estimators = names(all_estimator_results)
  for (estimator in estimators) {
    results = all_estimator_results[[estimator]]
    
    # Get full estimator name for the table
    full_estimator_name = switch(estimator,
                                  "IW" = "Iterative Weights",
                                  "NW" = "Normal Weights",
                                  "KDE" = "KDE",
                                  "FD" = "FD",
                                  estimator)
    
    cat(sprintf("%s & %.4f & %.4f & [%.4f, %.4f] \\\\\n", 
                full_estimator_name, 
                results$overall_slope, 
                results$theoretical_rate,
                results$conf_int[1], 
                results$conf_int[2]))
  }
  
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat(sprintf("\\label{tab:convergence_weibull_%s_%s}\n", current_shape, current_scale))
  cat("\\end{table}\n\n")
}

weibull_params = unique(all_results_df[, c("shape", "scale")])

for (i in 1:nrow(weibull_params)) {
  current_shape = weibull_params$shape[i]
  current_scale = weibull_params$scale[i]
  
  # Filter results for this shape/scale combination
  df_subset = subset(all_results_df, shape == current_shape & scale == current_scale)
  
  # Create a list to store results for all estimators
  all_estimator_results = list()
  
  # Analyze each estimator
  estimators = c("IW", "KDE", "NW", "FD")
  ise_columns = c("iter_ise", "kde_ise", "direct_ise", "fd_ise")
  
  for (e in 1:length(estimators)) {
    estimator = estimators[e]
    ise_col = ise_columns[e]
    
    # Analyze convergence
    analysis_results = analyze_convergence(df_subset, n_values, estimator, ise_col)
    all_estimator_results[[estimator]] = analysis_results
  }
  
  # Generate the LaTeX table for this parameter combination
  generate_latex_table_weibull(all_estimator_results, current_shape, current_scale)
}



saveRDS(results_all, file = "results_all.rds", compress = TRUE)
a =  readRDS("results_all.rds")
a
