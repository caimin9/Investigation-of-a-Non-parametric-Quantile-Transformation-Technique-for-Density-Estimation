library(dplyr)
library(MASS)
library(parallel)    # For detectCores()
library(doParallel)  # For registerDoParallel()
library(foreach)     # For foreach loops
options(pillar.sigfig = 4)

num_cores = detectCores() - 1  # leave one core free
cl = makeCluster(num_cores)
registerDoParallel(cl)

# Make sure worker processes can see necessary packages:
clusterEvalQ(cl, {
  library(dplyr)
  library(MASS)
})


run_combined_simulation_t = function(n, seed, dfParam, initial_df) {
  set.seed(seed)
  
  # 1) Generate t-dist data
  Y = rt(n, df = dfParam)
  
  # set up a symmetrical grid around 0
  upper_bound = qt(0.9999, df = dfParam)
  y_grid      = seq(-upper_bound, upper_bound, length.out = 1000)
  
  # "True" t-density
  true_density = dt(y_grid, df = dfParam)
  
  # 2) Common QDE setup
  probs    = ppoints(n)
  normal_q = qnorm(probs)
  
  spline_weights = function(quantiles, n, density_values) {
    quantiles * (1 - quantiles) / (n * (density_values^2))
  }
  
  trapz = function(x, y) {
    sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(x)]) / 2)
  }
  
  transformation = function(spline_obj, x_vals, f_x_in) {
    dy_dx = predict(spline_obj, x = x_vals, deriv = 1)$y
    f_x_in / abs(dy_dx)
  }
  
  #For the fitdistr
  my_t_density = function(x, m, s, df) {
    # Compute the density of the t-distribution with a location and scale adjustment.
    dt((x - m) / s, df = df) / s
  }
  #-------------------------------------------------------------
  # A) Iterative QDE
  #-------------------------------------------------------------
  sample_q = quantile(Y, probs)
  f_x      = dnorm(normal_q)
  
  iterative_qde = function(x_vals, y_vals, x_dist, probs, n, df_val, 
                           true_density, num_iter = 30) {
    weights        = spline_weights(probs, n, x_dist)
    weights_matrix = matrix(NA, nrow = n, ncol = num_iter + 1)
    error_vec      = numeric(num_iter)
    
    weights_matrix[,1] = weights
    
    for(iter in seq_len(num_iter)) {
      spline_fit  = smooth.spline(x_vals, y_vals, w = 1/weights, df = df_val)
      f_y_new     = transformation(spline_fit, x_vals, x_dist)
      weights     = spline_weights(probs, n, f_y_new)
      weights_matrix[,iter+1] = weights
      
      # Evaluate ISE at iteration
      y_pred_iter     = spline_fit$y
      f_y_interp_iter = approx(y_pred_iter, f_y_new, xout = y_grid, rule=2)$y
      err_sq          = (f_y_interp_iter - true_density)^2
      error_vec[iter] = trapz(y_grid, err_sq)
    }
    list(spline_fit     = spline_fit,
         final_f_y      = f_y_new,
         weights_matrix = weights_matrix,
         error_per_iter = error_vec)
  }
  
  # run iterative QDE with the user-chosen initial_df
  res_iter       = iterative_qde(normal_q, sample_q, f_x, probs, n, df_val = initial_df, 
                                 true_density, num_iter=30)
  best_iter_idx  = which.min(res_iter$error_per_iter)
  final_weights  = res_iter$weights_matrix[, best_iter_idx]
  
  # fine-tune the final df
  calc_ISE_iter = function(df_spline) {
    spline_fit = smooth.spline(normal_q, sample_q, w = 1/final_weights, df = df_spline)
    f_y_pred   = transformation(spline_fit, normal_q, f_x)
    f_y_interp = approx(spline_fit$y, f_y_pred, xout=y_grid, rule=2)$y
    trapz(y_grid, (f_y_interp - true_density)^2)
  }
  optim_iter      = optim(par=initial_df, fn=calc_ISE_iter, method="Brent", lower=1, upper=30)
  df_iter_opt     = optim_iter$par
  
  # final iterative estimate
  spline_iter_fin = smooth.spline(normal_q, sample_q, w = 1/final_weights, df = df_iter_opt)
  f_y_iter_fin    = transformation(spline_iter_fin, normal_q, f_x)
  iter_pdf        = approx(spline_iter_fin$y, f_y_iter_fin, xout=y_grid, rule=2)$y
  iter_ise        = trapz(y_grid, (iter_pdf - true_density)^2)
  
  #-------------------------------------------------------------
  # B) Direct QDE
  #-------------------------------------------------------------
  weights_direct = spline_weights(probs, n, f_x)
  
  calc_ISE_direct = function(df_spline) {
    spline_fit = smooth.spline(normal_q, sample_q, w = 1/weights_direct, df = df_spline)
    f_y_pred   = transformation(spline_fit, normal_q, f_x)
    f_y_interp = approx(spline_fit$y, f_y_pred, xout=y_grid, rule=2)$y
    trapz(y_grid, (f_y_interp - true_density)^2)
  }
  optim_direct     = optim(par=initial_df, fn=calc_ISE_direct, method="Brent", lower=1, upper=30)
  df_direct_opt    = optim_direct$par
  
  spline_direct_fin = smooth.spline(normal_q, sample_q, w = 1/weights_direct, df = df_direct_opt)
  f_y_direct_fin    = transformation(spline_direct_fin, normal_q, f_x)
  direct_pdf        = approx(spline_direct_fin$y, f_y_direct_fin, xout=y_grid, rule=2)$y
  direct_ise        = trapz(y_grid, (direct_pdf - true_density)^2)
  
  #-------------------------------------------------------------
  # C) KDE
  #-------------------------------------------------------------
  calc_ISE_kde = function(bw_val) {
    kde_fit = density(Y, kernel="gaussian", bw=bw_val,n=n)
    kde_pdf = approx(kde_fit$x, kde_fit$y, xout=y_grid, rule=2)$y
    trapz(y_grid, (kde_pdf - true_density)^2)
  }
  optim_kde   = optim(par=1, fn=calc_ISE_kde, method="Brent", lower=0.01, upper=4)
  bw_opt      = optim_kde$par
  
  kde_fit = density(Y, bw = bw_opt, kernel="gaussian", n=n)
  kde_pdf = approx(kde_fit$x, kde_fit$y, xout=y_grid, rule=2)$y
  kde_ise = trapz(y_grid, (kde_pdf - true_density)^2)
  
  #-------------------------------------------------------------
  # D) t-dist MLE (fitdistr)
  #-------------------------------------------------------------
  fit_norm = fitdistr(Y, densfun = 't',
                      lower = c(-Inf, 1e-06, 1),
                      upper = c(Inf,100, 4000)) 
  df_hat     = fit_norm$estimate["df"]
  
  fd_pdf     = dt(y_grid, df = df_hat)
  fd_pdf[is.infinite(fd_pdf)] = 1e-30
  fd_ise     = trapz(y_grid, (fd_pdf - true_density)^2)
  
  # Return everything
  list(
    n             = n,
    seed          = seed,
    dfParam       = dfParam,
    initial_df    = initial_df,
    iter_ise      = iter_ise,
    direct_ise    = direct_ise,
    kde_ise       = kde_ise,
    fd_ise        = fd_ise,   # normal MLE approach
    df_iter_opt   = df_iter_opt,
    df_direct_opt = df_direct_opt,
    bw_opt        = bw_opt,
    
    densities = list(
      y_grid       = y_grid,
      true_density = true_density,
      iter_density = iter_pdf,
      direct_density = direct_pdf,
      kde_density   = kde_pdf,
      fd_density    = fd_pdf
    )
  )
}


# Example sets of sample size, seeds, t-parameters (df), and initial df
n_values = 100 * 2^(0:7)
seeds = c(1:500)
t_params           = c(2.01,3,6, 12,20)   # degrees of freedom

initial_df_values  = c(4)        # initial spline df to try

# Master results container
#results_all = list()

for(dfp in t_params) {
  for(init_df in initial_df_values) {
    cat(sprintf("\nRunning parallel sim for t-dist with df=%.2f, initial_df=%g\n",
                dfp, init_df))
    
    # Create data frame for final results + a list for densities
    results_combined  = data.frame()
    density_examples  = list()
    
    # All combos of n, seed
    param_grid = expand.grid(n = n_values, seed = seeds)
    
    # 3.1) Parallel foreach loop
    parallel_results = foreach(i = 1:nrow(param_grid), .combine = rbind,
                                .packages = c("dplyr", "MASS")) %dopar% {
                                  n    = param_grid$n[i]
                                  seed = param_grid$seed[i]
                                  
                                  # optional progress message
                                  cat(sprintf("  Worker %d: n=%d, seed=%d\n", Sys.getpid(), n, seed))
                                  
                                  sim_res = run_combined_simulation_t(n = n, seed = seed,
                                                                       dfParam = dfp, initial_df = init_df)
                                  
                                  data.frame(
                                    dfParam     = dfp,
                                    initial_df  = init_df,
                                    n           = sim_res$n,
                                    seed        = sim_res$seed,
                                    iter_ise    = sim_res$iter_ise,
                                    direct_ise  = sim_res$direct_ise,
                                    kde_ise     = sim_res$kde_ise,
                                    fd_ise      = sim_res$fd_ise,  # normal MLE
                                    df_iter     = sim_res$df_iter_opt,
                                    df_direct   = sim_res$df_direct_opt,
                                    bw          = sim_res$bw_opt,
                                    # whether to store densities
                                    save_density= (seed == 1)
                                  )
                                }
    
    # 3.2) Combine parallel results
    results_combined = parallel_results
    
    # 3.3) Extract density examples for seed=1
    for(n in n_values) {
      sim_res = run_combined_simulation_t(
        n=n, seed=1, dfParam=dfp, initial_df=init_df
      )
      density_examples[[as.character(n)]] = sim_res$densities
    }
    
    # store in results_all
    key = sprintf("t_df_%.2f_initdf_%g", dfp, init_df)
    results_all[[key]] = list(
      results   = results_combined,
      densities = density_examples
    )
  }
}

stopCluster(cl)
# Suppose we pick one scenario: df=2.0, init_df=3
results_all = readRDS("t_results_all_merged.rds") #Do this if data already saved

chosen_key = "t_df_2.00_initdf_3"

if(chosen_key %in% names(results_all)) {
  results_chosen   = results_all[[chosen_key]]$results
  density_examples = results_all[[chosen_key]]$densities
  
  # 5.1) Boxplot of ISE by sample size
  n_unique = sort(unique(results_chosen$n))
  par(mfrow = c(1,2))
  for(nv in n_unique) {
    sub_df = results_chosen[results_chosen$n == nv, ]
    boxplot(sub_df$iter_ise, sub_df$direct_ise, sub_df$kde_ise, sub_df$fd_ise,
            names=c("Iter","Direct","KDE","MLE(norm)"),
            main = paste("n =", nv), ylab = "ISE")
  }
  
  # 5.2) Compare densities for largest n
  largest_n = max(n_unique)
  dens      = density_examples[[as.character(largest_n)]]
  
  plot(dens$y_grid, dens$true_density, type="l", col="black", lwd=2,
       main = sprintf("t-dist df=%.2f, init_df=%g, n=%d", 
                      results_chosen$dfParam[1], 
                      results_chosen$initial_df[1],
                      largest_n),
       xlab = "y", ylab = "Density")
  lines(dens$y_grid, dens$iter_density,   col="blue",   lty=2, lwd=2)
  lines(dens$y_grid, dens$direct_density, col="red",    lty=3, lwd=2)
  lines(dens$y_grid, dens$kde_density,    col="green",  lty=4, lwd=2)
  lines(dens$y_grid, dens$fd_density,     col="purple", lty=5, lwd=2)
  
  legend("topright",
         legend=c("True t-dist","IterQDE","DirectQDE","KDE","Normal MLE"),
         col   =c("black","blue","red","green","purple"),
         lty   =1:5, lwd=2, bty="n")
}

# 5.3) Summaries across all scenarios
all_results_df = do.call(rbind, lapply(results_all, function(x) x$results))


overall_summary = all_results_df %>%
  group_by(dfParam,n) %>%
  summarise(
    mean_iter = mean(iter_ise),
    mean_direct = mean(direct_ise),
    mean_kde = mean(kde_ise),
    mean_fd  = mean(fd_ise),
    sd_iter_ise   = sd(iter_ise),
    sd_direct_ise = sd(direct_ise),
    sd_kde_ise    = sd(kde_ise),
    sd_fd_ise     = sd(fd_ise),
    .groups  = "drop"
  )
print(overall_summary, n=1000)


names(all_results_df)


##########################################################################################
######### BOXPLOT OF ERRORS (SHOWCASE ERROR CHARCTERISTICS)
##########################################################################################


# Create a new plot showing log(ISE) vs log(n) with boxplots at each sample size point
par(mar = c(5, 5, 4, 2))  # Increase margins for better spacing

plot(NULL, 
     xlim = range(log(n_values)),
     ylim = range(log(results_combined$iter_ise)),
     main = 'Log of Integrated Squared Error vs Sample Size',
     xlab = 'log(Sample Size)',
     ylab = 'log(ISE)',
     cex.main = 1.5,    # Larger title
     cex.lab = 1.3,     # Larger axis labels
     cex.axis = 1.2,
     axes = FALSE)    # Larger axis numbers

# Now draw the x-axis with log(n) positions, but label them with n_values
axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)

# Draw the y-axis normally
axis(2)

# Add a box around the plot
box()



# Add boxplots at each sample size
boxplot_data = list()
for (i in 1:length(n_values)) {
  n_val = n_values[i]
  subset_data = results_combined[results_combined$n == n_val, ]
  log_ise_values = log(subset_data$iter_ise)
  
  # Add boxplot at the x position corresponding to log(n)
  boxplot(log_ise_values, add = TRUE, at = log(n_val), 
          boxwex = .7,     # Increase boxplot width
          col = "lightblue", 
          outline = TRUE, 
          outcol = "gray60", 
          outpch = 20, 
          outcex = 1,      # Larger outlier points
          axes = FALSE)
  
  # Store mean for the line
  boxplot_data[[i]] = list(n = n_val, mean_log_ise = mean(log_ise_values))
}

# Extract n values and mean log(ISE) values for the line
n_vals = sapply(boxplot_data, function(x) x$n)
mean_log_ise_vals = sapply(boxplot_data, function(x) x$mean_log_ise)

# Add a line connecting the means of the boxplots - make it thicker
lines(log(n_vals), mean_log_ise_vals, col = 'red', lwd = 3)
points(log(n_vals), mean_log_ise_vals, col = 'red', pch = 19, cex = 1.5)

# Add grid lines for better readability
#grid(lty = "dotted", col = "lightgray")

# Add a reference line showing n^(-1/2) convergence rate (optional)
n_ref = n_values[1]
ise_ref = exp(mean_log_ise_vals[1])
theoretical_line = log(ise_ref) - 4/5 * (log(n_values) - log(n_ref))
lines(log(n_values), theoretical_line, col = 'darkgreen', lwd = 2, lty = 2)

# Add legend with larger text and symbols
legend("topright", 
       legend = c("Mean log(ISE)", "Boxplot at each sample size", "n^(-4/5) rate"),
       col = c("red", "lightblue", "darkgreen"), 
       lwd = c(3, NA, 2), 
       lty = c(1, NA, 2),
       pch = c(19, NA, NA), 
       fill = c(NA, "lightblue", NA),
       bty = 'n',
       cex = 1.2)   # Larger legend text





###################################################################

# Create separate plots for each df parameter
#pdf("t_distribution_ise_plots.pdf", width = 10, height = 8)  # Optional: save to PDF

# Loop through each df parameter
for (df_param in t_params) {
  # Filter data for current df parameter
  df_subset = all_results_df[all_results_df$dfParam == df_param, ]
  
  # Calculate mean log(ISE) for each sample size
  mean_log_ise = df_subset %>%
    group_by(n) %>%
    summarise(mean_log_ise = mean(log(iter_ise)), .groups = "drop")
  
  # Sort by n to ensure proper line connection
  mean_log_ise = mean_log_ise[order(mean_log_ise$n), ]
  
  # Create plot
  par(mar = c(5, 5, 4, 2))  # Increase margins
  
  plot(log(mean_log_ise$n), mean_log_ise$mean_log_ise,
       type = "o", pch = 19, col = "blue", lwd = 2,
       main = sprintf("Log of ISE vs Sample Size - t-distribution (df = %.2f)", df_param),
       xlab = "log(Sample Size)", ylab = "log(ISE)",
       cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1,
       axes = FALSE)
  
  # Draw the axes with n_values
  axis(1, at = log(n_values), labels = n_values, cex.axis = 0.9)
  axis(2)
  box()
  
  # Add a reference line showing n^(-4/5) convergence rate
  n_ref = n_values[1]
  ise_ref = exp(mean_log_ise$mean_log_ise[mean_log_ise$n == n_ref])
  theoretical_line = log(ise_ref) - 10/11* (log(n_values) - log(n_ref))
  lines(log(n_values), theoretical_line, col = 'red', lwd = 2, lty = 2)
  
  # Add legend
  legend("topright", 
         legend = c(sprintf("df = %.2f", df_param), "n^(-4/5) rate"),
         col = c("blue", "red"), 
         lwd = c(2, 2), 
         lty = c(1, 2),
         pch = c(19, NA), 
         bty = 'n',
         cex = 1.1)
  
  # Optional: add boxplots at each sample size
  for (n_val in n_values) {
    n_subset = df_subset[df_subset$n == n_val, ]
    boxplot(log(n_subset$iter_ise), add = TRUE, at = log(n_val),
            boxwex = 0.15, col = "lightblue", outline = TRUE,
            outcol = "gray60", outpch = 20, outcex = 0.8, axes = FALSE)
  }
}



#############################################################################
###### ESTIMATOR ANALYSIS
#############################################################################

################################################################################
######## Relative Efficiency
##############################################################################

summary_table = all_results_df %>%
  group_by(dfParam,n) %>%
  summarize(
    mean_iter_vs_direct = mean(iter_ise/direct_ise),
    mean_direct_vs_kde = mean(direct_ise/kde_ise),
    mean_iter_vs_kde = mean(iter_ise/kde_ise)

  ) %>%
  ungroup()
print(summary_table, n = 100)



# Keep your existing data structure but correct the efficiency calculations
for(df_val in unique(all_results_df$dfParam)) {
  
  # Filter data for this specific distribution
  dist_data = all_results_df %>% 
    filter(dfParam == df_val)
  
  # Calculate mean relative efficiency correctly - first mean ISEs, then ratios
  mean_data = dist_data %>%
    group_by(n) %>%
    summarize(
      mean_iter_vs_direct = mean(direct_ise)/mean(iter_ise),
      mean_direct_vs_kde = mean(kde_ise)/mean(direct_ise),
      mean_iter_vs_kde = mean(kde_ise)/mean(iter_ise)
    )
  
  # Determine a good y-axis range with some padding
  y_min = min(c(mean_data$mean_iter_vs_direct, 
                 mean_data$mean_direct_vs_kde, 
                 mean_data$mean_iter_vs_kde)) * 0.9
  y_max = max(c(mean_data$mean_iter_vs_direct, 
                 mean_data$mean_direct_vs_kde, 
                 mean_data$mean_iter_vs_kde)) * 1.05
  
  # Set up improved plot parameters
  par(mar = c(5, 5, 4, 2), mgp = c(3, 1, 0), cex.main = 1.2)
  
  # Create proper t-distribution title with mathematical notation
  main_title = substitute(paste("Relative Efficiency for ", t(df)), 
                           list(df = df_val))
  
  # Create a line plot with all three metrics
  plot(mean_data$n, mean_data$mean_iter_vs_direct,
       type = "n", 
       log = "x",  # Log scale for x-axis
       xlim = range(mean_data$n),
       ylim = c(y_min, y_max),  # Dynamically adjusted y-axis
       xlab = "Sample Size", 
       ylab = "Mean Relative Efficiency",
       main = main_title,
       xaxt = "n") # Suppress default x-axis
  
  # Create custom x-axis with actual values
  axis(1, at = unique(mean_data$n), labels = unique(mean_data$n))
  
  # Add subtle grid for better readability
  grid(nx = NA, ny = NULL, lty = "dotted", col = "gray90")
  
  # Add reference line at y = 1
  abline(h = 1, lty = 2, col = "gray50")
  
  # Add the three lines with improved styling
  lines(mean_data$n, mean_data$mean_iter_vs_direct, 
        type = "o", col = "blue", pch = 16, lwd = 1.5)
  lines(mean_data$n, mean_data$mean_direct_vs_kde, 
        type = "o", col = "darkgreen", pch = 17, lwd = 1.5)
  lines(mean_data$n, mean_data$mean_iter_vs_kde, 
        type = "o", col = "red", pch = 15, lwd = 1.5)
  
  # Determine best legend position based on data patterns
  n_points = length(mean_data$n)
  last_few = max(1, n_points - 2):n_points  # Last 3 points or all if fewer
  
  # Check if upper right is crowded (high values at end of series)
  upper_right_values = c(
    mean_data$mean_iter_vs_direct[last_few],
    mean_data$mean_direct_vs_kde[last_few],
    mean_data$mean_iter_vs_kde[last_few]
  )
  
  # Check if upper left is crowded (high values at start of series)
  upper_left_values = c(
    mean_data$mean_iter_vs_direct[1:min(3, n_points)],
    mean_data$mean_direct_vs_kde[1:min(3, n_points)],
    mean_data$mean_iter_vs_kde[1:min(3, n_points)]
  )
  
  # Determine legend position based on data density
  if(max(upper_right_values, na.rm = TRUE) > 0.8 * y_max) {
    # Upper right is crowded with high values
    if(max(upper_left_values, na.rm = TRUE) > 0.8 * y_max) {
      # Both corners are crowded, use bottom right
      legend_pos = "bottomright"
    } else {
      # Upper left has space
      legend_pos = "topleft"
    }
  } else {
    # Upper right has space
    legend_pos = "topright"
  }
  
  # Add transparent legend with correct labels
  legend(legend_pos, 
         legend = c("NW vs Iterative", "KDE vs NW", "KDE vs Iterative"),
         col = c("blue", "darkgreen", "red"), 
         pch = c(16, 17, 15), 
         lty = 1,
         lwd = 1.5,
         pt.cex = 1.1,
         y.intersp = 1.2,  # Increased spacing between items
         cex = 1.0,
         bty = "n")        # Transparent
}

#######################################################################
########## KDE vs Iterative for t-distributions
###########################################################################
comparison_data = all_results_df %>%
  group_by(dfParam, n) %>%
  summarize(
    mean_kde_vs_iter = mean(kde_ise)/mean(iter_ise),  # Corrected calculation
    .groups = "drop"
  )

# Set up improved plot parameters
par(mar = c(5, 5, 4, 2), mgp = c(3, 1, 0), cex.main = 1.2)

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_kde_vs_iter),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = expression("KDE vs Iterative Relative Efficiency"),
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add subtle grid for better readability
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray90")

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray50")

# Use Zissou 1 color palette
dist_colors = hcl.colors(length(unique(comparison_data$dfParam)), "Zissou 1")
dist_counter = 1

# Add lines for each distribution
for(df_val in unique(comparison_data$dfParam)) {
  dist_data = comparison_data %>% 
    filter(dfParam == df_val) %>%
    arrange(n) # Make sure data is ordered by sample size
  
  lines(dist_data$n, dist_data$mean_kde_vs_iter, 
        type = "o", 
        col = dist_colors[dist_counter], 
        pch = 15 + dist_counter,
        lwd = 1.5)  # Slightly thicker lines
  
  dist_counter = dist_counter + 1
}

# Build labels with t-distribution notation
legend_labels = paste0("t(", unique(comparison_data$dfParam), ")")

# Determine best legend position based on data patterns
y_values_end = tapply(comparison_data$mean_kde_vs_iter[comparison_data$n == max(comparison_data$n)], 
                       comparison_data$dfParam[comparison_data$n == max(comparison_data$n)], mean)

if(mean(y_values_end, na.rm=TRUE) > 0.8 * max(comparison_data$mean_kde_vs_iter)) {
  legend_pos = "topleft"
} else {
  legend_pos = "topright"
}

legend("topleft", 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       lwd = 1.5,
       pt.cex = 1.1,
       y.intersp = 1.2,
       bty = "n")

#######################################################################
########## KDE vs NW for t-distributions
###########################################################################
comparison_data = all_results_df %>%
  group_by(dfParam, n) %>%
  summarize(
    mean_kde_vs_nw = mean(kde_ise)/mean(direct_ise),  # Corrected calculation
    .groups = "drop"
  )

# Set up improved plot parameters
par(mar = c(5, 5, 4, 2), mgp = c(3, 1, 0), cex.main = 1.2)

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_kde_vs_nw),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = expression("KDE vs Normal Weights Relative Efficiency"),
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add subtle grid for better readability
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray90")

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray50")

# Use Zissou 1 color palette
dist_colors = hcl.colors(length(unique(comparison_data$dfParam)), "Zissou 1")
dist_counter = 1

# Add lines for each distribution
for(df_val in unique(comparison_data$dfParam)) {
  dist_data = comparison_data %>% 
    filter(dfParam == df_val) %>%
    arrange(n) # Make sure data is ordered by sample size
  
  lines(dist_data$n, dist_data$mean_kde_vs_nw, 
        type = "o", 
        col = dist_colors[dist_counter], 
        pch = 15 + dist_counter,
        lwd = 1.5)  # Slightly thicker lines
  
  dist_counter = dist_counter + 1
}

# Build labels with t-distribution notation
legend_labels = paste0("t(", unique(comparison_data$dfParam), ")")

# Determine best legend position based on data patterns
y_values_end = tapply(comparison_data$mean_kde_vs_nw[comparison_data$n == max(comparison_data$n)], 
                       comparison_data$dfParam[comparison_data$n == max(comparison_data$n)], mean)

if(mean(y_values_end, na.rm=TRUE) > 0.8 * max(comparison_data$mean_kde_vs_nw)) {
  legend_pos = "topleft"
} else {
  legend_pos = "topright"
}

legend('topleft', 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       lwd = 1.5,
       pt.cex = 1.1,
       y.intersp = 1.2,
       bty = "n")

#######################################################################
########## Iter vs Direct
###########################################################################
#######################################################################
########## NW vs Iterative for t-distributions
###########################################################################
comparison_data = all_results_df %>%
  group_by(dfParam, n) %>%
  summarize(
    mean_nw_vs_iter = mean(direct_ise)/mean(iter_ise),  # Corrected calculation
    .groups = "drop"
  )

# Set up improved plot parameters
par(mar = c(5, 5, 4, 2), mgp = c(3, 1, 0), cex.main = 1.2)

# Create an empty plot with proper x-axis values
plot(NULL, 
     xlim = range(comparison_data$n), 
     ylim = range(comparison_data$mean_nw_vs_iter),
     log = "x",
     xlab = "Sample Size", 
     ylab = "Mean Relative Efficiency",
     main = expression("Normal Weights vs Iterative Relative Efficiency"),
     xaxt = "n") # Suppress default x-axis

# Create custom x-axis with actual values
axis(1, at = unique(comparison_data$n), labels = unique(comparison_data$n))

# Add subtle grid for better readability
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray90")

# Add reference line at y = 1
abline(h = 1, lty = 2, col = "gray50")

# Use Zissou 1 color palette
dist_colors = hcl.colors(length(unique(comparison_data$dfParam)), "Zissou 1")
dist_counter = 1

# Add lines for each distribution
for(df_val in unique(comparison_data$dfParam)) {
  dist_data = comparison_data %>% 
    filter(dfParam == df_val) %>%
    arrange(n) # Make sure data is ordered by sample size
  
  lines(dist_data$n, dist_data$mean_nw_vs_iter, 
        type = "o", 
        col = dist_colors[dist_counter], 
        pch = 15 + dist_counter,
        lwd = 1.5)  # Slightly thicker lines
  
  dist_counter = dist_counter + 1
}

# Build labels with t-distribution notation
legend_labels = paste0("t(", unique(comparison_data$dfParam), ")")

# Determine best legend position based on data patterns
y_values_end = tapply(comparison_data$mean_nw_vs_iter[comparison_data$n == max(comparison_data$n)], 
                       comparison_data$dfParam[comparison_data$n == max(comparison_data$n)], mean)

if(mean(y_values_end, na.rm=TRUE) > 0.8 * max(comparison_data$mean_nw_vs_iter)) {
  legend_pos = "topleft"
} else {
  legend_pos = "topright"
}

legend('bottomleft', 
       legend = legend_labels,
       col = dist_colors, 
       pch = 16:(15 + length(legend_labels)), 
       lty = 1,
       lwd = 1.5,
       pt.cex = 1.1,
       y.intersp = 1.2,
       bty = "n")




##############################################################################################################
##############                        ESTIMATION CHARACTERISTICS
##############################################################################################################

############################################################################################################
############### Log of iter ISE VS n with Larger Boxplots
############################################################################################################

par(mfrow = c(1,1))
t_params = unique(all_results_df$dfParam)
for (df_param in t_params) {
  
  df_subset = subset(all_results_df, dfParam == df_param)
  
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
       main = substitute(paste("Log of Iterative ISE vs Sample Size - ", "t(", s, ")"), 
                         list(s = df_param)),
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
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 1  # Position above the line
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
         legend = c(bquote(t(.(df_param))), 
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
t_params = unique(all_results_df$dfParam)
for (df_param in t_params) {
  
  df_subset = subset(all_results_df, dfParam == df_param)
  
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
       main = substitute(paste("Log of KDE ISE vs Sample Size - ", "t(", s, ")"), 
                         list(s = df_param)),
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
         legend = c(bquote(t(.(df_param))), 
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
t_params = unique(all_results_df$dfParam)
for (df_param in t_params) {
  
  df_subset = subset(all_results_df, dfParam == df_param)
  
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
  
  # Same colour scheme as before
  box_color = "#B3D9FF"  # Light blue for boxplots
  line_color = "#0066CC"  # Darker blue for the line
  point_color = "#003366"  # Even darker blue for points
  ref_color = "#CC3366"    # Pink/red for reference line
  slope_color = "#228B22"  # Green for slope values
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$direct_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = substitute(paste("Log of NW ISE vs Sample Size - ", "t(", s, ")"), 
                         list(s = df_param)),
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
    mid_y = (y_values[j] + y_values[j+1]) / 2 + 1 # Increased vertical offset
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
         legend = c(bquote(t(.(df_param))), 
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
t_params = unique(all_results_df$dfParam)
for (df_param in t_params) {
  
  df_subset = subset(all_results_df, dfParam == df_param)
  
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
  
  # Same colour scheme as before
  box_color = "#B3D9FF"  # Light blue for boxplots
  line_color = "#0066CC"  # Darker blue for the line
  point_color = "#003366"  # Even darker blue for points
  ref_color = "#CC3366"    # Pink/red for reference line
  slope_color = "#228B22"  # Green for slope values
  
  # Improve the layout and margins
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Set up the plot with appropriate axis limits
  y_range = range(log(df_subset$fd_ise), na.rm = TRUE)
  y_padding = 0.1 * diff(y_range)
  
  # Create empty plot with grid
  plot(log(summary_stats$n), summary_stats$mean_log_ise,
       type = "n",
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = substitute(paste("Log of FD ISE vs Sample Size - ", "t(", s, ")"), 
                         list(s = df_param)),
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
  
  # Create a better looking text box for the slope info - MODIFIED SIZE
  # Keep at bottom center but make wider
  text_box_x = (par("usr")[1] + par("usr")[2]) / 2  # Center of x-axis
  text_box_y = par("usr")[3] + .8 * diff(par("usr")[3:4]) * 0.15  # Just above bottom margin
  
  # Make the box much wider to fit all text
  box_width = 4.5  # Increase this value to make the box wider
  box_height = 1  # Increase this value to make the box taller
  
  # Create a wider text box
  rect(text_box_x - box_width,  
       text_box_y - box_height/2, 
       text_box_x + box_width, 
       text_box_y + box_height/2, 
       col = rgb(1, 1, 1, 0.8), 
       border = "gray70")
  
  # Add the slope value text with the same position
  text(text_box_x, text_box_y, 
       sprintf("Overall slope: %.3f (Theoretical: -0.800)", overall_slope),
       cex = 0.9, col = "black")
  
  # Add legend with no box (bty='n')
  legend("topright", 
         legend = c(bquote(t(.(df_param))), 
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
##############################################################################################################
# Function to analyze convergence rates
analyze_convergence = function(df_subset, n_values, estimator_name = "IW", ise_column = "iter_ise") {
  # Calculate mean for log(ISE)
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(!!sym(ise_column))),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  # Calculate the slope using linear regression
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit linear model to get overall slope
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Return analysis results
  return(list(
    summary_stats = summary_stats,
    x_values = x_values,
    y_values = y_values,
    overall_slope = overall_slope
  ))
}

# Function to create comparison plots with all estimators
create_comparison_plot = function(df_subset, n_values, df_param) {
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
  
  # Get y-range for all estimators
  all_y_values = unlist(lapply(all_estimator_results, function(x) x$y_values))
  y_range = range(all_y_values)
  y_padding = 0.1 * diff(y_range)
  
  # Set up clean plot with grid
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Create plot area
  plot(log(n_values), log(n_values), type = "n",
       xlim = range(log(n_values)),
       ylim = c(min(y_range) - y_padding, max(y_range) + y_padding),
       main = paste0("Comparison of Estimators - t(", df_param, ")"),
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
  
  # Colors and line types for the estimators
  colors = c("#0066CC",    # Blue for Iterative Weights
              "#009933",    # Green for KDE
              "#FF6600",    # Orange for Normal Weights
              "#CC3366")    # Purple for FD (MLE)
  
  # Line types (solid, dashed, dotted, etc.)
  ltys = c(1, 2, 3, 4)
  
  # Point types (circle, triangle, square, etc.)
  pch_types = c(16, 17, 15, 18)
  
  # Plot each estimator's line
  for (e in 1:length(estimators)) {
    estimator = estimators[e]
    results = all_estimator_results[[estimator]]
    
    # Add line
    lines(results$x_values, results$y_values, 
          col = colors[e], lwd = 2, lty = ltys[e])
    
    # Add points
    points(results$x_values, results$y_values,
           pch = pch_types[e], col = colors[e], cex = 1.2)
  }
  
  # Add legend with full estimator names and slopes
  legend_text = sapply(1:length(estimators), function(e) {
    estimator = estimators[e]
    results = all_estimator_results[[estimator]]
    
    full_name = switch(estimator,
                        "IW" = "Iterative Weights",
                        "KDE" = "KDE",
                        "NW" = "Normal Weights",
                        "FD" = "FD",
                        estimator)
    
    sprintf("%s (slope: %.3f)", full_name, results$overall_slope)
  })
  
  # Add legend in top right
  legend("topright", 
         legend = legend_text,
         col = colors, 
         lwd = 2, 
         lty = ltys, 
         pch = pch_types,
         pt.cex = 1.2, 
         bty = 'n', 
         cex = 1.0)
}

# Create comparison plots for each t-distribution parameter
t_params = unique(all_results_df$dfParam)

for (df_param in t_params) {
  # Get subset of results for current parameter
  df_subset = subset(all_results_df, dfParam == df_param)
  
  # Create the comparison plot with all estimators
  create_comparison_plot(df_subset, n_values, df_param)
}









#############################################################
##### LAtex Table
##########################################################

analyze_convergence = function(df_subset, n_values, estimator_name = "IW", ise_column = "iter_ise") {
  # Calculate mean log(ISE) by sample size
  summary_stats = df_subset %>%
    group_by(n) %>%
    summarise(
      mean_log_ise = mean(log(!!sym(ise_column))),
      .groups = "drop"
    ) %>%
    arrange(n)
  
  x_values = log(summary_stats$n)
  y_values = summary_stats$mean_log_ise
  
  # Fit a linear model to obtain the empirical rate (slope)
  lm_fit = lm(y_values ~ x_values)
  overall_slope = coef(lm_fit)[2]
  
  # Extract the 95% confidence interval for the slope (second coefficient)
  conf_int = confint(lm_fit)[2, ]
  
  # The theoretical convergence rate for n^(-4/5) is -0.8 (if that is what you expect)
  theoretical_rate = -0.8
  
  # Return all results, including new fields
  return(list(
    summary_stats = summary_stats,
    x_values = x_values,
    y_values = y_values,
    overall_slope = overall_slope,
    conf_int = conf_int,
    theoretical_rate = theoretical_rate
  ))
}
# Function to generate a LaTeX table of results for the t-distribution
generate_latex_table_t = function(all_estimator_results, df_param) {
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{Convergence Rate Analysis for t(", df_param, ")}\n")
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
  cat(sprintf("\\label{tab:convergence_t_%s}\n", df_param))
  cat("\\end{table}\n\n")
}

# Generate LaTeX tables for each t-distribution parameter
for (df_param in t_params) {
  # Filter results for this df parameter
  df_subset = subset(all_results_df, dfParam == df_param)
  
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
  generate_latex_table_t(all_estimator_results, df_param)
}



##########
# Save results
saveRDS(results_all, file = "t_results_all.rds", compress = TRUE)
