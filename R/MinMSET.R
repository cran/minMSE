# Authors
# Sebastian Schneider, sschneider@coll.mpg.de; sebastian@sebastianschneider.eu
# Giulia Baldini, giulia.baldini@uni-bonn.de

# Copyright (C) 2019 Sebastian Schneider & Giulia Baldini

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


############################################################################
#### Min MSE Treatment Assignment for One or Multiple Treatment Groups #####
############################################################################
library(MASS)

plotting_file <- "plotting_data.csv"
output_file <- "output_file.txt"

evaluate_solution_matrix <- function(treatment, data, mse_weights) {
  ind <- c(1:length(treatment))
  n_treatments <- max(treatment)
  n_outcomes <- nrow(mse_weights)
  v_norm <- sum(apply(mse_weights, 2, vector_gcd))
  w_norm <- sum(apply(mse_weights, 1, vector_gcd))

  x_0 <- data[(treatment[ind] == 0), , drop = FALSE]
  x_0_inv <- ginv(crossprod(x_0)) # Import MASS and use ginv for the pseudoinverse
  inv_x0 <- w_norm * v_norm * x_0_inv # Weigh x_0 by the w and v norm

  wk_weight <- 1
  vk_weight <- 1
  sum_of_scaled <- inv_x0 # Initialize sum with inv_x0
  for (i in 1:n_outcomes){
    sum_of_inv_xt <- 0
    wk_weight <- vector_gcd(mse_weights[i, ])
    for (j in 1:n_treatments) {
      x <- data[(treatment[ind] == j), , drop = FALSE]
      x_inv <- ginv(crossprod(x))
      vk_weight <- mse_weights[i, j] / wk_weight
      sum_of_inv_xt <- sum_of_inv_xt + vk_weight * x_inv
    }
    sum_of_scaled <- sum_of_scaled + wk_weight * sum_of_inv_xt
  }

  mean_of_xi <- colMeans(data) # Mean of all columns
  evaluated_value <- crossprod(mean_of_xi, sum_of_scaled) %*% mean_of_xi
  return (evaluated_value)
}


evaluate_solution_vector <- function(treatment, data, mse_weights) {
  ind <- c(1:length(treatment))
  n_treatments <- max(treatment)
  v_norm <- sum(mse_weights) # Compute the l1 norm of the vector
  
  x_0 <- data[(treatment[ind] == 0), , drop = FALSE]
  x_0_inv <- ginv(crossprod(x_0)) # Import MASS and use ginv for the pseudoinverse
  inv_x0 <- v_norm * x_0_inv # Weigh x_0 by the v norm
  
  vk_weight <- 1
  sum_of_inv_xt <- inv_x0 # Initialize sum with inv_x0
  for (j in 1:n_treatments) {
    x <- data[(treatment[ind] == j), , drop = FALSE]
    x_inv <- ginv(crossprod(x))
    vk_weight <- mse_weights[j]
    sum_of_inv_xt <- sum_of_inv_xt + vk_weight * x_inv
  }
  
  mean_of_xi <- colMeans(data) # Mean of all columns
  evaluated_value <- crossprod(mean_of_xi, sum_of_inv_xt) %*% mean_of_xi
  
  return (evaluated_value)
}

evaluate_solution <- function(treatment, data, mse_weights = NULL) {
  ind <- c(1:length(treatment))
  n_treatments <- max(treatment)

  x_0 <- data[(treatment[ind] == 0), , drop = FALSE]
  x_0_inv <- ginv(crossprod(x_0)) # Import MASS and use ginv for the pseudoinverse
  sum_of_inv_xt <- n_treatments * x_0_inv # Weigh x_0 by the number of treatments
  
  for (i in 1:n_treatments) {
    x <- data[(treatment[ind] == i), , drop = FALSE]
    x_inv <- ginv(crossprod(x))
    sum_of_inv_xt <- sum_of_inv_xt + x_inv
  }
  
  mean_of_xi <- colMeans(data) # Mean of all columns
  evaluated_value <- crossprod(mean_of_xi, sum_of_inv_xt) %*% mean_of_xi

  return (evaluated_value)
}


evaluate_solution.optim <- function(par, data, evaluation_function = evaluate_solution, swap_treatment_function = NULL, mse_weights = NULL, change = NULL, prev_index_list = NULL) { # Optim wants fn and gr to have the same parameters
  return(evaluation_function(par, data, mse_weights))
}

scale_vars <- function(data) {
  col_means <- colMeans(data, na.rm = TRUE)
  data_sd <- sapply(data, sd, na.rm = TRUE)
  for (i in 1:ncol(data)) {
    data[is.na(data[, i]), i] <- col_means[i]
    if (data_sd[i] > 0) {
      data[, i] <- scale(data[, i], center = FALSE, scale = data_sd[i])
    }
    else {
      data[, i] <- rep(0 , length(data[, i]))
    }
  }
  return (data)
}

swap_treatment <- function(current_treatment, change, prev_index_list = NULL) {
  max_swaps <- change - 1 # Maximum possible swaps between indices
  s_ind <- sample(1:length(current_treatment)) # Sample random indices so that we can always swap the first max_swaps elements
  scrambled_treatment <- current_treatment[s_ind]
  scrambled_treatment[1:(max_swaps+1)] <- scrambled_treatment[c(2:(max_swaps+1), 1)]
  current_treatment <- scrambled_treatment[order(s_ind)] # If we don't have a previous assignment, take the whole list
  
  return (current_treatment)
}

swap_treatment_prev <- function(current_treatment, change, prev_index_list) {
  max_swaps <- min(change, length(prev_index_list)) - 1 # Maximum swaps is either change or the length of the ones still to assign - it might be that there are less to assign than 'change'
  trunc_treatment <- current_treatment[prev_index_list] # Treatment containing only the ones to assign
  s_ind <- sample(1:length(trunc_treatment)) # Sample random indices so that we can always swap the first max_swaps elements
  scrambled_treatment <- trunc_treatment[s_ind]
  scrambled_treatment[1:(max_swaps+1)] <- scrambled_treatment[c(2:(max_swaps+1), 1)]
  current_treatment[prev_index_list] <- scrambled_treatment[order(s_ind)] # Otherwise substitute only the ones that are not assigned
  
  return (current_treatment)
}

swap_treatment.optim <- function(current_treatment, data = NULL, evaluation_function = NULL, swap_treatment_function = swap_treatment, mse_weights = NULL, change, prev_index_list = NULL) { # Optim wants fn and gr to have the same parameters
  return (swap_treatment_function(current_treatment, change, prev_index_list))
}

sample_with_prev_treatment <- function(prev_treatment, n_treatments, n_per_group){
  groups = n_treatments + 1;
  to_assign <- sum(is.na(prev_treatment)) # How many people still have to be assigned
  current_treatment <- prev_treatment
  remaining_times <- c()
  for (i in 1:groups){
    if(is.na(n_per_group[2])){ # The number of elements per group is numerical
      current_n_per_group <- n_per_group
    } else { # The number of elements per group is a vector
      current_n_per_group <- n_per_group[i]
    }
    remaining_times[i] <- current_n_per_group - length(which(current_treatment == (i - 1))) # vector that contains how many times each index can still be assigned
  }
  smaller_repetitions <- rep(0:n_treatments, times = remaining_times)
  sampled <- sample(smaller_repetitions, size = to_assign, replace = FALSE)

  j <- 1

  for (i in 1:length(current_treatment)){
    if(is.na(current_treatment[i])){
      current_treatment[i] <- sampled[j]
      j <- j + 1
    }
  }
  return (current_treatment)
}

assign_treatment <- function(current_data,
                             prev_treatment = NULL,
                             evaluation_function = evaluate_solution,
                             swap_treatment_function = swap_treatment,
                             n_treatments = 1,
                             n_per_group = NULL,
                             mse_weights = NULL,
                             iterations = 50,
                             change = 3,
                             cooling = 1,
                             t0 = 10,
                             tmax = 10,
                             built_in = 0,
                             plot = 0,
                             create_plot_file = 1) {

  n_obs <- dim(current_data)[1]
  if(is.null(n_per_group)){
    n_per_group <- ceiling(n_obs / (n_treatments + 1))
    repetitions <- rep(0:n_treatments, each = n_per_group)
  } else { # in case division is user provided
    repetitions <- rep(0:n_treatments, n_per_group)
  }

  # Generate first solution
  if (is.null(prev_treatment)){
    current_treatment <- sample(repetitions, size = nrow(current_data), replace = FALSE)
    prev_treat_na_index <- NULL
  } else { # In case some people are already assigned a group
    current_treatment <- sample_with_prev_treatment(prev_treatment, n_treatments, n_per_group)
    prev_treat_na_index <- which(is.na(prev_treatment)) # Contains the indices of the non-assigned people
  }

  # Evaluate first solution
  current_ssd <- evaluation_function(current_treatment, current_data, mse_weights)

  # Choose the best solution out of 5% of the iterations
  num_first_solutions = min(max(round(iterations / 100 * 5) - 1, n_obs - 1), round(iterations / 100 * 10) - 1)
  
  for (k in 1:num_first_solutions) {
    if (is.null(prev_treatment)){
      next_treatment <- sample(repetitions, size = nrow(current_data), replace = FALSE)
    } else { # In case some people are already assigned a group
      next_treatment <- sample_with_prev_treatment(prev_treatment, n_treatments, n_per_group)
    }
    next_ssd <- evaluation_function(next_treatment, current_data, mse_weights)

    if (next_ssd < current_ssd) {
      current_treatment <- next_treatment
      current_ssd <- next_ssd
    }
  }

  if (t0 < 0) {
    t0 <- -(1 / t0) * current_ssd
  }

  opt_treatment <- current_treatment
  opt_ssd <- current_ssd
  if (plot){
    opts <- c()
  }

  if (plot && create_plot_file) {
    if (built_in){
      vec <- 0:iterations
      vec <- vec[seq(1, length(vec), 10)]
      vec <- append(vec, iterations - 1, after = length(vec) - 1)
      plotting_data <- data.frame(iter = vec)
    } else {
      plotting_data <- data.frame(iter = 0:iterations)
    }
    write.csv(plotting_data, paste(tempdir(), plotting_file, sep="/"))
  }

  # Use the built-in function instead of our program
  if (built_in == 1) {
    save_output_filename <- paste(tempdir(), output_file, sep="/")
    sink(save_output_filename)
    number <- as.double(regmatches(opt_ssd, regexec('-?[0-9\\.]+', opt_ssd))[[1]])
    custom_scale <- opt_ssd / number
    result <-
      optim(
        par = current_treatment,
        fn = evaluate_solution.optim,
        evaluation_function = evaluation_function,
        swap_treatment_function = swap_treatment_function,
        mse_weights = mse_weights,
        data = current_data,
        gr = swap_treatment.optim,
        change = change,
        prev_index_list = prev_treat_na_index,
        method = 'SANN',
        control = c(
          fnscale = custom_scale, # Depends on the data, it allows printing
          trace = 1,
          maxit = iterations,
          temp = t0,
          tmax = tmax,
          REPORT = 1
        )
      )

    opt_treatment <- result$par
    opt_ssd <- result$value

    if (built_in && plot) {
      out_file <- readLines(save_output_filename)
      k <- 1
      for (line in out_file) {
        y <- regmatches(line, regexec('value (-?[0-9\\.]+)', line))
        if (length(y[[1]]) != 0) {
          opts[k] <- as.double(strsplit(y[[1]], " ")[[2]]) * custom_scale
          k <- k + 1
        }
      }
      plotting_data <- read.csv(paste(tempdir(), plotting_file, sep="/"))
      plotting_data[ , ncol(plotting_data) + 1] <- opts
      plotting_data <- plotting_data[,-1]
      write.csv(plotting_data, paste(tempdir(), plotting_file, sep="/"))
    }

    sink()

    return (append(opt_treatment, opt_ssd))
  }

  if (plot) {
    opts[1] <- opt_ssd
  }

  # Start optimization with our program
  for (k in 1:iterations) {
    # Swap some number (change) of values in the treatment assignment
    next_treatment <- swap_treatment_function(current_treatment, change, prev_treat_na_index)
    next_ssd <- evaluation_function(next_treatment, current_data, mse_weights)

    rand = runif(1)
    if (cooling == 1) {
      temperature <- t0 / log(floor((k - 1) / tmax) * tmax + exp(1))
    } else {
      temperature <- t0 / (floor((k - 1) / tmax) * tmax + 1)
    }

    # Simulated Annealing
    if (rand <=  min(1, exp((current_ssd - next_ssd) / (temperature)))) {
      current_treatment <- next_treatment
      current_ssd <- next_ssd

      if (current_ssd < opt_ssd) {
        opt_treatment <- current_treatment
        opt_ssd <- current_ssd
      }
    }

    if (plot) {
      opts[k + 1] <- opt_ssd
    }
  }

  if (plot) {
    plotting_data <- read.csv(paste(tempdir(), plotting_file, sep="/"))
    plotting_data[ , ncol(plotting_data) + 1] <- opts
    plotting_data <- plotting_data[,-1]
    write.csv(plotting_data, paste(tempdir(), plotting_file, sep="/"))
  }

  return (append(opt_treatment, opt_ssd))
}

assign_minMSE_treatment <- function(data,
                                    prev_treatment = NULL,
                                    n_treatments = 1,
                                    n_per_group = NULL,
                                    mse_weights = NULL,
                                    iterations = 50,
                                    change = 3,
                                    cooling = 1,
                                    t0 = 10,
                                    tmax = 10,
                                    built_in = 0,
                                    desired_test_vectors = 100,
                                    percentage_equal_treatments = 1,
                                    plot = 0,
                                    trace_output = 1,
                                    filename = NULL) {
  if (is.matrix(data)){ # If the data is a matrix
   data <- as.data.frame(data) # convert to dataframe
  } else if (!is.data.frame(data)){
    stop("Please convert the data to a dataframe or to a matrix before calling this function.")
  }

  n_obs <- dim(data)[1]

  # Checking all variables
  if (n_treatments <= 0) {
    stop("Number of treatments must be greater or equal than 1.\n")
  }

  if (!is.null(n_per_group)){
    if ((is.vector(n_per_group) && length(n_per_group) != (n_treatments+1)) || !is.vector(n_per_group)) {
      stop("Group sizes must be provided in a vector with the same length as the number of experimental groups!\n")
    }
    
    if (sum(n_per_group) != n_obs){
      stop("The group sizes have to be consistent with the total number of observations.") 
    }
  }
  

  if (!is.null(mse_weights) && !is.vector(mse_weights) && !is.matrix(mse_weights)){
    stop("The weights must be either a vector or a matrix.")
  }

  if (!is.null(mse_weights) &&
      ((is.vector(mse_weights) && n_treatments != length(mse_weights)) ||
       (is.matrix(mse_weights) && n_treatments != ncol(mse_weights)))){
    stop("The number of treatment mse_weights is different from the number of treatments.")
  }

  if (iterations < 1) {
    stop("Positive number of iterations needed.\n")
  }

  if (change < 1) {
    change <- 3
  }

  if (change >= n_obs) {
    change <- min(3, n_obs)
  }

  if (cooling != 2 && cooling != 1) {
    stop("The value of cooling should be either 1 (cooling function 1) or 2 (cooling function 2). If built_in is 1, then the cooling function is 1.\n")
  }

  if (t0 == 0) {
    t0 <- 10
  }

  if (tmax < 1) {
    tmax <- 10
  }

  if (built_in != 0 && built_in != 1) {
    stop("The value of built_in should be either 0 or 1 (use built_in R function).\n")
  }

  if (desired_test_vectors < 1) {
    stop("Positive maximum number of tests needed.\n")
  }

  if (plot != 0 && plot != 1) {
    stop("The value of plot should be either 0 (no plotting) or 1 (plotting).\n")
  }

  if (trace_output != 0 && trace_output != 1) {
    stop("The value of trace_output should be either 0 (no trace) or 1 (with trace).\n")
  }

  if (percentage_equal_treatments > 100 || percentage_equal_treatments < 1){
    stop("The percentage of equal treatments should be a value between 1 and 100.\n")
  }

  if (!is.null(prev_treatment) && length(prev_treatment) != nrow(data)){
    stop("The length of the previous treatment should be the same as the number of observations.")
  }
  
  if (is.vector(mse_weights)){
    evaluation_function <- evaluate_solution_vector
  } else if (is.matrix(mse_weights)){
    evaluation_function <- evaluate_solution_matrix
  } else {
    evaluation_function <- evaluate_solution
  }
  
  if (is.null(prev_treatment)){
    swap_treatment_function <- swap_treatment
  } else {
    swap_treatment_function <- swap_treatment_prev
  }

  # Set missing values to a variable's mean and rescale (without centering) all variables to have a standard deviation/variance of 1
  current_data = scale_vars(data)

  # Run the treatment assignment many times for testing
  curr_iter <- 0
  curr_max_equal <- 0
  df_treatments <- data.frame(1:n_obs)

  if (plot) {
    if (built_in){
      vec <- 0:iterations
      vec <- vec[seq(1, length(vec), 10)]
      vec <- append(vec, iterations - 1, after = length(vec) - 1)
      plotting_data <- data.frame(iter = vec)
    } else {
      plotting_data <- data.frame(iter = 0:iterations)
    }
    write.csv(plotting_data, paste(tempdir(), plotting_file, sep="/"))
  }

  current_count_equal <- 0 # We keep track of how many treatment vectors are the same
  max_equal_treatments <- (percentage_equal_treatments * desired_test_vectors) / 100 # Exact number of treatments that are allowed to be equal
  while (curr_iter < desired_test_vectors) {
    curr_iter <- curr_iter + 1
    curr_treatment <- assign_treatment(as.matrix(current_data), prev_treatment, evaluation_function, swap_treatment_function, n_treatments, n_per_group, mse_weights, iterations, change, cooling, t0, tmax, built_in, plot, 0)

    curr_opt_ssd <- tail(curr_treatment, n=1)
    curr_treatment <- curr_treatment[-length(curr_treatment)]
    if (percentage_equal_treatments < 100){
      current_count_equal <- current_count_equal + count_occurrences(df_treatments, curr_treatment) # If the new vector is equal to others, we increase the count
      if (current_count_equal > max_equal_treatments) {
        message ("Number of allowed duplicates exceeded; if you need more, consider decreasing 'iterations'.\n")
        curr_iter <- curr_iter - 1
        break
      }
    }

    iterstring <- paste('treatment_iter_', curr_iter, sep = '')

    df_treatments[, iterstring] <- curr_treatment

    # Print some feedback output
    if (trace_output) {
      cat ("Test iteration no. ", curr_iter, " of ", desired_test_vectors, ".\n", sep = "")
    }
  }

  if (desired_test_vectors == 1){
    message ("Value of the objective function for the currently minimal MSE treatment assignment found: ", curr_opt_ssd, "\n\n", sep = "")
  } else {
    message ('The maximum number of iterations with a confidence of ', percentage_equal_treatments, '% is: ', curr_iter, '.\n\n', sep = '')
  }

  if (plot) {
    plotting_data <- read.csv(paste(tempdir(), plotting_file, sep="/"))
    plotting_data <- plotting_data[, -1]

    dev.new()
    iter_vec <- plotting_data[, 'iter']
    curr_col = curr_iter + 1
    min_val <- min(plotting_data[, 2:curr_col])
    max_val <- max(plotting_data[, 2:curr_col])

    plot(0, 0, xlim = c(min(iter_vec), max(iter_vec)), ylim = c(min_val, max_val), type = "n",
         ylab = 'Value of objective function', xlab = 'Iterations', main = 'Optimization of the treatment assignment')
    cl <- rainbow(curr_iter)
    for (i in 1:curr_iter) {
      lines(iter_vec, plotting_data[, i + 1], col = cl[i])
    }
  }

  names(df_treatments)[1] <- "id"

  if (!is.null(filename)){
    if(grepl("/", filename)){
      write.csv(df_treatments, filename)
    }else {
      write.csv(df_treatments, paste(tempdir(), filename, sep="/"))
    }
  }

  return (df_treatments)
}
