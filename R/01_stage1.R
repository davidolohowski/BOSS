#' Bayesian Optimization via Mode-Seeking Surrogate (BOSS)
#'
#' Run a GP-based sequential design algorithm that (every \code{update_step} iter):
#' \enumerate{
#'   \item Re-estimates GP length-scale/noise;
#'   \item Optimizes an acquisition function (UCB) to pick the next \code{x};
#'   \item Optionally checks a modal‐based convergence criterion.
#' }
#'
#' @param func A function \code{f: [0,1]^D -> R} to optimize.
#' @param update_step Integer; how often to refit the GP hyperparameters.
#' @param max_iter Maximum number of BO iterations.
#' @param D Input dimension.
#' @param lower Numeric vector of length D giving lower bounds (in original scale).
#' @param upper Numeric vector of length D giving upper bounds.
#' @param nu GP covariance smoothness. \eqn{\nu = \infty} represents squared-exponential kernel;
#' \eqn{3 < \nu < \infty} represents the Matern kernel with smoothness \eqn{\lceil\nu\rceil - 1/2}.
#' \eqn{\nu > 3} is required to ensure convergence properties.
#' @param quad Logical; if \code{FALSE} (default), use a zero-mean GP model.
#'   If \code{TRUE}, fit a linear-plus-quadratic mean model
#'   \eqn{\beta_0 + x^\top \beta_1 + (x \otimes x)^\top \beta_2} alongside the GP.
#' @param noise_var GP observation noise variance.
#' @param modal_iter_check How often (in iterations) to check convergence via modal criteria.
#' @param modal_check_warmup Minimum iter before starting modal checks.
#' @param modal_k.nn Number of nearest neighbours to use in modal criterion.
#' @param modal_eps Convergence threshold for the modal criterion.
#' @param initial_design Number of initial points to sample.
#' @param delta Exploration parameter passed to UCB.
#' @param optim.n Number of multistart in acquisition optimization.
#' @param optim.max.iter Maxit for length-scale optimization.
#' @param opt.lengthscale.grid Optional grid size for length-scale search.
#' @param opt.grid Optional grid size for acquisition search.
#' @param verbose Verbosity level (0–3).
#'
#' @return A "boss_modal" object (S3 class \code{boss_modal}), which is a list containing:
#'   \describe{
#'     \item{\code{objective_function}}{The original objective function \code{f}.}
#'     \item{\code{D}}{Input dimension of the objective function.}
#'     \item{\code{design_points}}{List with components \code{x} (scaled inputs), \code{x_original} (original inputs), and \code{y} (responses).}
#'     \item{\code{gp_params}}{List with \code{length_scale} and \code{signal_var}, the final GP hyperparameters.}
#'     \item{\code{surrogate}}{A function \code{(xvalue) -> numeric} giving the GP posterior mean at \code{xvalue}.}
#'     \item{\code{essential_support}}{Placeholder for essential support (to be computed later).}
#'     \item{\code{fill_in}}{Placeholder for additional outputs (to be computed later).}
#'     \item{\code{modal_result}}{Data frame of modal-difference diagnostics collected during optimization.}
#'   }
#'
#' @importFrom lhs randomLHS
#' @importFrom optimx multistart
#' @importFrom numDeriv hessian
#' @importFrom MASS mvrnorm
#' @export
BOSS_modal <- function(func, update_step = 5, max_iter = 100, D = 1,
                       quad = FALSE,
                       lower = rep(0, D), upper = rep(1, D),
                       nu = Inf,
                       noise_var = 1e-6,
                       modal_iter_check = 10,  modal_check_warmup = 20,
                       modal_k.nn = 5, modal_eps = 0.1, # The number of nearest neighbor K should be smaller than modal_iter_check, otherwise there may not be enough changes in the mode
                       initial_design = 5, delta = 0.01,
                       optim.n = 5, optim.max.iter = 1000,
                       opt.lengthscale.grid = NULL,  # Grid-based option for lengthscale optimization.
                       opt.grid = NULL,              # Grid-based option for AF optimization.
                       verbose = 3) {              # Verbosity level: 0, 1, 2, or 3

  # Initialize a helper for verbose printing
  vprint <- function(level, msg) {
    if (verbose >= level) print(msg)
  }

  # If verbose >= 1, print that stage 1 starts.
  if (verbose >= 1) {
    cat("Stage 1: Bayesian Optimization via Mode-Seeking Surrogate (BOSS) started.\n")
    start_time <- Sys.time()
  }

  # Check if dimensions of lower and upper bounds match.
  if(length(lower) != D || length(upper) != D) {
    stop("lower and upper must have the same length as the function's input dimension")
  }

  if(nu <= 3){
    stop('GP kernel smoothness parameter must be nu > 3.')
  }

  # Initialize matrices/vectors to store evaluations.
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()

  if (verbose == 3) {
    print('Initial evaluation phase...')
  }

  if(D > 1){
    if (verbose == 3) {
      print("Using Latin Hypercube Sampling for initial design for D > 1.")
    }
    initial <- lhs::randomLHS(initial_design, D)
    initial <- t(apply(initial, 1, function(x) x*(upper - lower) + lower))
  }else{
    if (verbose == 3) {
      print("Using equally spaced spaced initial design for D = 1.")
    }
    initial <- matrix(rep(seq(from = (lower), to = (upper), length.out = initial_design), D),
                      nrow = initial_design, ncol = D, byrow = FALSE)
  }

  for (i in 1:nrow(initial)) {
    xmat_trans <- rbind(xmat_trans, (initial[i,] - lower)/(upper - lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }
  # Center the function values.
  rel <- mean(yvec)
  yvec <- yvec - rel
  num_initial <- initial_design

  # Initial optimization for lengthscale.
  if (!is.null(opt.lengthscale.grid)) {
    length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
    like_vec <- sapply(length_scale_vec, function(l)
      -compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   quad = quad, D = D, noise_var = noise_var, nu = nu))
    max_idx <- which.max(like_vec)
    length_scale <- length_scale_vec[max_idx]
    lik <- like_vec[max_idx]
    signal_var <- compute_like(length_scale = NULL, y = yvec, x = xmat_trans, quad = quad, D = D)
  } else {
    opt <- optim(runif(D, 0.01, 0.99), function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   quad = quad, D = D, noise_var = noise_var, nu = nu),
      control = list(maxit = optim.max.iter), lower = rep(0.01, D), upper = rep(0.99, D), method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
    signal_var <- compute_like(length_scale = NULL, y = yvec, x = xmat_trans, quad = quad, D = D)
  }
  vprint(3, paste("The new length.scale:", length_scale))
  vprint(3, paste("The new signal_var:", signal_var))

  choice_cov <- cov_generator(length_scale = length_scale, signal_var = signal_var, nu = nu)

  # Initialize modal-related quantities.
  modal_max_diff <- Inf
  modal_f_old <- list(mode = NA, hessian = NA)
  modal_result <- data.frame(modal_max_diff = numeric(0), i = numeric(0))

  # Initialize the grid for the acquisition function, if provided.
  if (!is.null(opt.grid)) {
    grid_list <- replicate(D, seq(0, 1, length.out = opt.grid), simplify = FALSE)
    grid <- as.matrix(expand.grid(grid_list))
  }

  i <- 1
  while(i <= max_iter && modal_eps < modal_max_diff){
    if(verbose == 3) {
      print(paste("Iteration:", i))
    } else if(verbose == 2) {
      cat(paste("Iteration:", i, "\n"))
    }

    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

    if(i %% update_step == 0){
      if(verbose == 3) print("Time to update the parameters!")
      if(!is.null(opt.lengthscale.grid)) {
        length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
        like_vec <- sapply(length_scale_vec, function(l)
          -compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                        quad = quad, D = D, noise_var = noise_var, nu = nu))
        max_idx <- which.max(like_vec)
        length_scale <- length_scale_vec[max_idx]
        lik <- like_vec[max_idx]
        signal_var <- compute_like(length_scale = NULL, y = newdata$y, x = newdata$x, quad = quad, D = D)
      } else {
        opt <- optim(runif(D, 0.01, 0.99), function(l)
          compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                       quad = quad, D = D, noise_var = noise_var, nu = nu),
          control = list(maxit = optim.max.iter), lower = rep(0.01, D), upper = rep(0.99, D), method = 'L-BFGS-B')
        length_scale <- opt$par
        lik <- opt$value
        signal_var <- compute_like(length_scale = NULL, y = newdata$y, x = newdata$x, quad = quad, D = D)
      }
      if(verbose == 3) {
        print(paste("The new length.scale:", length_scale))
        print(paste("The new signal_var:", signal_var))
      }
      choice_cov <- cov_generator(length_scale = length_scale, signal_var = signal_var, nu = nu)
    }

    # Optimize the acquisition function (UCB).

    # Use pre-computed quantities for the acquisition function.
    intern <- predict_gp_internal(data = newdata,
                                  noise_var = noise_var,
                                  choice_cov = choice_cov,
                                  quad = quad)

    # UCB_internal(x, intern, covfn, nv, D, d)

    if(verbose == 3) print("Maximize Acquisition Function")
    if(is.null(opt.grid)){
      # Multi-start local optimization.
      initialize_UCB <- matrix(runif(D*optim.n, rep(0, D), rep((1), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB_internal(x = matrix(x, nrow = 1, ncol = D), intern = intern, D = D,
                     covfn = choice_cov, nv = noise_var, num_training = i + num_initial, delta = delta),
        control = list(maxit = 100),
        lower = rep(0, D), upper = rep((1), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      af_values <- apply(grid, 1, function(x)
        -UCB_internal(x = matrix(x, nrow = 1, ncol = D), intern = intern, D = D,
                      covfn = choice_cov, nv = noise_var, num_training = i + num_initial, delta = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # Take out the best point from the grid.
      grid <- grid[-best_idx, , drop = FALSE]
    }
    # if next_point is already covered, skip
    if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
      if(verbose == 3){
        # ("Next point is already covered, skipping.")
        print(paste("Next point:", next_point*(upper - lower) + lower, "is already covered, skipping."))
      }
      next
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)

    yvec <- c(yvec, y_new - rel)

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with modal every modal_iter_check iterations.
    if(i %% modal_iter_check == 0  && i >= modal_check_warmup){
      if(verbose == 3) print("Time to check modal difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var, quad = quad)$mean
      }

      lf_design <- list(x = xmat_trans, y = yvec)
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lf_design))
      # find current optimizer
      mode_point <- xmat_trans[which.max(yvec),]

      # find hessian at the current optimizer
      mode_hess <- numDeriv::hessian(fn_new, matrix(mode_point, nrow = 1, ncol = D))
      mode_hess <- (mode_hess + t(mode_hess))/2

      # compute the trace of the hessian
      post_sd <- 1/sqrt(sum(abs(eigen(mode_hess)$values)))

      # Compute Euclidean distances from the best point to all other design points.
      distances <- sqrt(rowSums((xmat_trans - matrix(rep(mode_point, nrow(xmat_trans)),
                                                     nrow = nrow(xmat_trans), byrow = TRUE))^2))
      # Exclude the zero distance (self-distance)
      nonzero_distances <- distances[distances > 0]

      # find average nearest neighbor distance
      nn_dists <- mean(nonzero_distances[order(nonzero_distances)[1:modal_k.nn]])

      #stop if nn_dists is less than certain percentage (acc) of post_sd
      mode <- nn_dists/post_sd
      modal_f_new <- list()
      modal_f_new$mode <- mode
      modal_f_new$hessian <- 1/post_sd

      modal_diff_new <- list(
        mode = modal_f_new$mode,
        hessian = abs(modal_f_new$hessian - modal_f_old$hessian)/abs(modal_f_old$hessian)
      )

      # if any of the modal differences is NaN or NA, set it to Inf
      modal_diff_new$mode <- ifelse(is.na(modal_diff_new$mode) | is.nan(modal_diff_new$mode), Inf, modal_diff_new$mode)
      modal_diff_new$hessian <- ifelse(is.na(modal_diff_new$hessian) | is.nan(modal_diff_new$hessian), Inf, modal_diff_new$hessian)

      if(verbose == 3) {
        print(paste("Modal rel-difference:", modal_diff_new$mode))
        print(paste("Hessian rel-difference in second moment:", modal_diff_new$hessian))
      }

      modal_max_diff <- max(modal_diff_new$mode, modal_diff_new$hessian)
      modal_diff <- modal_diff_new
      modal_f_old <- modal_f_new
      modal_result <- rbind(modal_result, data.frame(modal_max_diff = modal_max_diff, i = i))
    }
    i <- i + 1
  }

  if (modal_max_diff < modal_eps && modal_eps > 0) {
    if(verbose >= 2) print("Posterior surrogate converged based on modal criteria!")
  } else {
    if(verbose >= 2) {
      print(paste0("Maximum iterations reached! Maximum modal difference: ", modal_max_diff))
    }
  }

  # If verbose >= 1, print that stage 1 ends.
  if (verbose >= 1) {
    cat("Stage 1: BOSS finished.\n")
    end_time <- Sys.time()
    cat("Total time taken:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
  }

  intern <- predict_gp_internal(
    data = list(x = xmat_trans, y = (yvec)),
    noise_var = noise_var,
    choice_cov = choice_cov,
    quad = quad
  )

  surrogate_fn <- function(xvalue) {
    xvalue <- matrix(xvalue, ncol = D)
    xvalue_transform <- sweep(xvalue, 2, lower, "-") / (upper - lower)

    obtain_mean_internal(
      xnew = xvalue_transform, intern = intern,
      covfn = choice_cov,
      D = D) + rel
  }

  boss_result <- structure(list(
    objective_function = func, D = D,
    lower = lower, upper = upper, quad = quad, noise_var = noise_var,
    design_points = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
    gp_params = list(length_scale = length_scale, signal_var = signal_var, nu = nu),
    surrogate = surrogate_fn,
    essential_support = NULL,  # This will be computed in the next stage.
    fill_in = NULL, # This will be computed in the next stage.
    modal_result = modal_result),
    class = "boss_modal")

  return(boss_result)
}


#' Bayesian Optimization via Sequential MCMC (In Development)
#'
#' Run a GP-based sequential MCMC algorithm that:
#' \enumerate{
#'   \item Re-estimates GP length-scale/noise;
#'   \item Alternate between performing BO (\code{inner_iter} steps) and perform MCMC sampling on the current GP surrogate to update BO search region;
#'   \item Algorithm terminates when successive MCMC samples indicate convergence (by \eqn{\hat{R}} statistics) and BO has reached modal convergence.
#' }
#'
#' @param func A function \code{f: [0,1]^D -> R} to optimize.
#' @param update_step Integer; how often to refit the GP hyperparameters.
#' @param max_iter Maximum number of BO iterations.
#' @param inner_iter Integer; how often to perform MCMC sampling based on the current surrogate
#' @param MCMC_size Integer; MCMC iterations
#' @param D Input dimension.
#' @param lower Numeric vector of length D giving lower bounds (in original scale).
#' @param upper Numeric vector of length D giving upper bounds.
#' @param nu GP covariance smoothness. \eqn{\nu = \infty} represents squared-exponential kernel;
#' \eqn{3 < \nu < \infty} represents the Matern kernel with smoothness \eqn{\lceil\nu\rceil - 1/2}.
#' \eqn{\nu > 3} is required to ensure convergence properties.
#' @param quad Logical; if \code{FALSE} (default), use a zero-mean GP model.
#'   If \code{TRUE}, fit a linear-plus-quadratic mean model
#'   \eqn{\beta_0 + x^\top \beta_1 + (x \otimes x)^\top \beta_2} alongside the GP.
#' @param noise_var GP observation noise variance.
#' @param modal_k.nn Number of nearest neighbours to compute modal covering.
#' @param initial_design Number of initial points to sample.
#' @param delta Exploration parameter passed to UCB.
#' @param alpha Numeric value between \code{(0,1)}. Used to determine new search region based on current MCMC sample.
#' @param acc Numeric; accuracy required of the modal exploitation by BO.
#' @param explore_size Integer; how many MCMC sample to use to explore location with the highest surrogate uncertainty.
#' @param optim.n Number of multistart in acquisition optimization.
#' @param optim.max.iter Maxit for length-scale optimization.
#' @param Rhat_eps Numeric; required accuracy to terminate the algorithm based on the \eqn{\hat{R}} statistics between successive MCMC samples.
#' Default to \code{Rhat_eps = 0.01} to represent a stopping criterion of \eqn{\hat{R} < 1.01}.
#' @param UCB_prob Numeric; probability threshold to indicate the modal convergence from UCB. Default to \code{UCB_prob < 0.1}.
#' @param verbose Verbosity level (0–3).
#'
#' @return A "boss_mcmc" object (S3 class \code{boss_mcmc}), which is a list containing:
#'   \describe{
#'     \item{\code{objective_function}}{The original objective function \code{f}.}
#'     \item{\code{D}}{Input dimension of the objective function.}
#'     \item{\code{design_points}}{List with components \code{x} (scaled inputs), \code{x_original} (original inputs), and \code{y} (responses).}
#'     \item{\code{gp_params}}{List with \code{length_scale} and \code{signal_var}, the final GP hyperparameters.}
#'     \item{\code{surrogate}}{A function \code{(xvalue) -> numeric} giving the GP posterior mean at \code{xvalue}.}
#'     \item{\code{essential_support}}{A box specifying the essential support identified by the last MCMC sample.}
#'     \item{\code{fill_in}}{Placeholder for additional outputs (to be computed later).}
#'     \item{\code{mcmc_result}}{List with components \code{sample} (MCMC sample from the last iteration), \code{R_hat} (\eqn{\hat{R}} statistics across MCMC-sampling phase.) }
#'     \item{\code{modal_result}}{Data frame of modal-distance diagnostics collected during optimization.}
#'   }
#'
#' @importFrom lhs randomLHS
#' @importFrom optimx multistart
#' @importFrom numDeriv hessian
#' @importFrom MASS mvrnorm
#' @importFrom coda mcmc.list
#' @importFrom coda as.mcmc
#' @importFrom coda gelman.diag
BOSS_mcmc <- function(func,
                      update_step = 10,
                      max_iter = 50,
                      inner_iter = 10,
                      MCMC_size = 10000,
                      D = 1,
                      lower = rep(0,D), upper = rep(1,D),
                      nu = Inf,
                      quad = FALSE,
                      noise_var = 1e-6,
                      modal_k.nn = 5,
                      initial_design = 10,
                      delta = 0.01,
                      alpha = 0.05, acc = 0.25,
                      explore_size = 500,
                      optim.n = 1, optim.max.iter = 1000,
                      Rhat_eps = 0.01, UCB_prob = 0.1,
                      verbose = 1) {

  # Initialize a helper for verbose printing
  vprint <- function(level, msg) {
    if (verbose >= level) print(msg)
  }

  # If verbose >= 1, print that stage 1 starts.
  if (verbose >= 1) {
    cat("Stage 1: Bayesian Optimization via Sequential MCMC started.\n")
    start_time <- Sys.time()
  }

  # Check if dimensions of lower and upper bounds match.
  if(length(lower) != D || length(upper) != D) {
    stop("lower and upper must have the same length as the function's input dimension")
  }

  if(nu <= 3){
    stop('GP kernel smoothness parameter must be nu > 3.')
  }

  eps <- Rhat_eps
  prob <- UCB_prob

  # Store the original bounds (used for normalization in the GP)
  original_lower <- lower
  original_upper <- upper

  ## INITIAL DESIGN: Sample from the original region and normalize using original bounds
  xmat <- NULL          # Evaluations in original scale
  xmat_trans <- NULL    # Normalized evaluations (relative to original bounds)
  yvec <- c()

  if(verbose == 3){
    cat("Initial evaluation phase...\n")
  }

  initial <- lhs::randomLHS(initial_design, D)
  initial <- t(apply(initial, 1, function(x) x*(original_upper - original_lower) + original_lower))

  for (i in 1:nrow(initial)) {
    xmat_trans <- rbind(xmat_trans, (initial[i,] - original_lower) / (original_upper - original_lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }

  # Center the observations
  rel <- mean(yvec)
  yvec <- yvec - rel

  if (verbose == 3) {
    cat("  Updating GP hyperparameters...\n")
  }

  opt <- optim(runif(D, 0.01, 0.99),
               function(l) compute_like(length_scale = l, x = xmat_trans, y = yvec,
                                        noise_var = noise_var,
                                        D, quad = quad, nu = nu),
               control = list(maxit = 100),
               lower = rep(0.01, D), upper = rep(0.99, D), method = "L-BFGS-B")
  current_length_scale <- opt$par
  current_signal_var <- compute_like(length_scale = NULL, y = yvec, x = xmat_trans,
                                     noise_var = noise_var,
                                     D, quad = quad, nu = nu)
  lik <- opt$value
  if (verbose == 3){
    cat("    GP log likelihood:", lik, "\n")
    cat("    New length_scale:", current_length_scale, "\n")
    cat("    New signal_var:", current_signal_var, "\n")
  }

  # The "search region" used for acquisition and SPG injection.
  # Initially, this is the same as the original region.
  choice_cov <- cov_generator(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)

  modal_result <- data.frame(modal_dist = numeric(0), UCB_prob = numeric(0), iter = numeric(0))
  mcmc_result <- data.frame(R_hat = numeric(0), iter = numeric(0))

  if (verbose == 2){
    cat("Initial MCMC run... \n")
  }
  newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

  intern_GP <- predict_gp_internal(
    data = newdata,
    noise_var = noise_var,
    choice_cov = choice_cov,
    quad = quad
  )


  surrogate <- function(xvalue) {
    sim_once_internal(
      xnew = xvalue, intern = intern_GP,
      covfn = choice_cov,
      D = D)
  }

  # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
  MCMC_sample <- adaptive_MCMC_internal(N_sample = MCMC_size, D = D, BO_surrogate = surrogate,
                                        start_point = newdata$x[which.max(newdata$y), ])

  # Convert the normalized MCMC samples back to the original scale using original bounds.
  UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
  UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower

  if (verbose == 3) {
    cat("  New search region updated.\n")
    cat("  New lower bound:", UCB_lower, "\n")
    cat("  New upper bound:", UCB_upper, "\n")
  }

  area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
  if (verbose == 3) {
    cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")
  }

  # Update the search region used for candidate generation
  lower_current <- UCB_lower
  upper_current <- UCB_upper

  old_MCMC_sample <- MCMC_sample
  R_hat <- Inf

  UCB_counter <- 0

  ## MAIN BO LOOP (total BO iterations = max_iter)
  for (iter in 1:max_iter) {
    if(verbose == 3) {
      print(paste("BO Iteration:", i))
    } else if(verbose == 2) {
      cat(paste("BO Iteration:", i, "\n"))
    }

    # Update GP hyperparameters every update_step iterations (using the GP data normalized by original bounds)
    if (iter %% update_step == 0) {
      if (verbose == 3) {
        cat("  Updating GP hyperparameters...\n")
      }


      opt <- optim(runif(D, 0.01, 0.99),
                   function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
                                            noise_var = noise_var,
                                            D, quad = quad, nu = nu),
                   control = list(maxit = 100),
                   lower = rep(0.01, D), upper = rep(0.99, D), method = "L-BFGS-B")
      current_length_scale <- opt$par
      current_signal_var <- compute_like(length_scale = NULL, y = yvec, x = xmat_trans,
                                         noise_var = noise_var,
                                         D, quad = quad, nu = nu)
      lik <- opt$value
      if (verbose == 3) {
        cat("    GP log likelihood:", lik, "\n")
        cat("    New length_scale:", current_length_scale, "\n")
        cat("    New signal_var:", current_signal_var, "\n")
      }
    }

    # Build the covariance function using the current hyperparameters
    choice_cov <- cov_generator(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)

    # Compute the adaptive probability for UCB based on the current design points.
    # Use the point with the best observed value (proxy for the mode).
    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

    intern_GP <- predict_gp_internal(
      data = newdata,
      noise_var = noise_var,
      choice_cov = choice_cov,
      quad = quad
    )


    surrogate <- function(xvalue) {
      obtain_mean_internal(
        xnew = xvalue, intern = intern_GP,
        covfn = choice_cov,
        D = D)
    }


    mode_point <- xmat_trans[which.max(yvec),]
    mode_hess <- numDeriv::hessian(surrogate, matrix(mode_point, nrow = 1, ncol = D))
    mode_hess <- (mode_hess + t(mode_hess))/2
    c_const <- sqrt(sum(abs(eigen(mode_hess)$values)))/acc*log(3)

    # Compute Euclidean distances from the best point to all other design points.
    distances <- sqrt(rowSums((xmat_trans - matrix(rep(mode_point, nrow(xmat_trans)),
                                                   nrow = nrow(xmat_trans), byrow = TRUE))^2))
    # Exclude the zero distance (self-distance)
    nonzero_distances <- distances[distances > 0]
    d_min <- mean(nonzero_distances[order(nonzero_distances)[1:modal_k.nn]])

    modal_distance <- round(d_min*sqrt(sum((original_upper - original_lower)^2)), 6)

    if (verbose == 3) {
      cat('  Current nearest neighbor distance around mode:', modal_distance, '\n')
    }

    p_value <- 2 / (1 + exp(-c_const * d_min)) - 1
    if (verbose == 3) {
      cat("  Adaptive UCB probability:", round(p_value, 4), "\n")
    }

    modal_result <- rbind(modal_result, data.frame(modal_dist = modal_distance, UCB_prob = p_value, iter = iter))

    UCB_flag <- runif(1) < p_value
    if (verbose == 3) {
      cat('  Do UCB?', UCB_flag, '\n')
    }
    if(UCB_flag){
      lower_current <- UCB_lower
      upper_current <- UCB_upper

      # Maximize the acquisition function (here UCB) to select the next evaluation point.
      # For acquisition, we sample candidates from the current search region.
      # When normalizing for the GP, we always use the original bounds.
      init_point <- lhs::randomLHS(optim.n, D)
      init_point <- t(apply(init_point, 1, function(x) x*(upper_current - lower_current) + lower_current))
      init_point_trans <- t(apply(init_point, 1, function(x) (x - original_lower)/(original_upper - original_lower)))

      if (verbose == 3) {
        cat("  Maximizing Acquisition Function...\n")
      }
      optimizer <- optimx::multistart(init_point_trans,
                                      function(x) UCB_internal(matrix(x, nrow = 1, ncol = D),
                                                               intern = intern_GP, covfn = choice_cov, nv = noise_var, D = D,
                                                               num_training = UCB_counter + 1, delta = delta),
                                      control = list(maxit = optim.max.iter),
                                      lower = (lower_current - original_lower) / (original_upper - original_lower),
                                      upper = (upper_current - original_lower) / (original_upper - original_lower),
                                      method = "L-BFGS-B")

      next_point_trans <- unlist(unname(unique(optimizer[which.min(optimizer$value), c(1:D)])))
    }
    else{
      thinned <- unique(MCMC_sample[-c(1:1000),])
      s_thinned <- thinned[sample(1:nrow(thinned), explore_size),]

      if (verbose == 3) {
        cat("  Maximizing Acquisition Function...\n")
      }
      var_thinned <- -apply(s_thinned, 1, function(x) EXPLORE_internal(matrix(x, nrow = 1, ncol = D),
                                                                       intern = intern_GP,
                                                                       covfn = choice_cov, nv = noise_var, D = D))

      next_point_trans <- s_thinned[which.max(var_thinned),]
    }

    UCB_counter <- UCB_counter + UCB_flag
    next_point <- next_point_trans * (original_upper - original_lower) + original_lower
    f_next <- func(next_point)

    # Append the new evaluation: for the GP we normalize using the original bounds.
    xmat_trans <- rbind(xmat_trans, (next_point - original_lower) / (original_upper - original_lower))
    xmat <- rbind(xmat, next_point)
    yvec <- c(yvec, f_next - rel)

    if (verbose == 3) {
      cat("  Next point:", paste(round(next_point, 4), collapse = ", "), "\n")
      cat("  Function value:", round(f_next, 4), "\n")
    }

    ## MCMC UPDATE: Every inner_iter BO iterations, update the search region via MCMC.
    if (iter %% inner_iter == 0) {
      if(verbose == 3) {
        cat("  Running MCMC update at iteration", iter)
      } else if(verbose == 2) {
        cat("  Running MCMC update at iteration", iter, '\n')
      }

      newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

      intern_GP <- predict_gp_internal(
        data = newdata,
        noise_var = noise_var,
        choice_cov = choice_cov,
        quad = quad
      )


      surrogate <- function(xvalue) {
        sim_once_internal(
          xnew = xvalue, intern = intern_GP,
          covfn = choice_cov,
          D = D)
      }

      # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
      MCMC_sample <- adaptive_MCMC_internal(N_sample = MCMC_size, D = D, BO_surrogate = surrogate,
                                            start_point = newdata$x[which.max(newdata$y), ])

      # Convert the normalized MCMC samples back to the original scale using original bounds.
      UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
      UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower


      if (verbose == 3) {
        cat("  New search region updated.\n")
        cat("  New lower bound:", UCB_lower, "\n")
        cat("  New upper bound:", UCB_upper, "\n")
      }
      area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
      if (verbose == 3) {
        cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")
      }

      # Compute convergence diagnostic (R_hat)
      chains <- coda::mcmc.list(coda::as.mcmc(MCMC_sample[-c(1:1000), ]),
                                coda::as.mcmc(old_MCMC_sample[-c(1:1000), ]))
      R_hat <- abs(max(coda::gelman.diag(chains)$psrf[, "Point est."]) - 1)

      mcmc_result <- rbind(mcmc_result, data.frame(R_hat = R_hat + 1, iter = iter))

      if (verbose == 3) {
        cat("  R_hat of MCMC samples:", R_hat + 1, "\n")
      }
      old_MCMC_sample <- MCMC_sample

      if (R_hat < eps & p_value < prob) {
        if (verbose >= 2) {
          cat("Convergence achieved with R_hat =", round(R_hat, 5), "and UCB probability =", round(p_value, 4), "\n")
        }
        break
      }
    }
  }

  if ((R_hat >= eps | p_value >= prob) & verbose >= 2) {
    cat("Maximum iterations reached and convergence criteria not met (R_hat =", round(R_hat, 5), "UCB probability =", round(p_value, 4), ")\n")
  }

  # If verbose >= 1, print that stage 1 ends.
  if (verbose >= 1) {
    cat("Stage 1: BOSS-MCMC finished.\n")
    end_time <- Sys.time()
    cat("Total time taken:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
  }

  newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

  intern_GP <- predict_gp_internal(
    data = newdata,
    noise_var = noise_var,
    choice_cov = choice_cov,
    quad = quad
  )

  surrogate_fn <- function(xvalue) {
    xvalue <- matrix(xvalue, ncol = D)
    xvalue_transform <- sweep(xvalue, 2, lower, "-") / (upper - lower)

    obtain_mean_internal(
      xnew = xvalue_transform, intern = intern_GP,
      covfn = choice_cov,
      D = D) + rel
  }

  MCMC_sample <- t(apply(MCMC_sample, 1, function(x) x*(original_upper - original_lower) + original_lower))

  boss_result <- structure(list(
    objective_function = func, D = D,
    lower = original_lower, upper = original_upper, quad = quad, noise_var = noise_var,
    design_points = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
    gp_params = list(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu),
    surrogate = surrogate_fn,
    essential_support = data.frame(lower = apply(MCMC_sample, 2, min),
                                   upper = apply(MCMC_sample, 2, max)),
    fill_in = NULL, # This will be computed in the next stage.
    mcmc_result = list(sample = MCMC_sample, R_hat = mcmc_result),
    modal_result = modal_result),
    class = "boss_mcmc")

  return(boss_result)
}


#' Update Mode Hessian for a \code{boss_modal} Object (with PSD check)
#'
#' Compute the numerical Hessian of the objective at the current best point
#' (in original input space), verify it is negative semi-definite (NSD), and
#' store it back into the \code{boss_result}.
#'
#' @param boss_result A S3 \code{boss_modal} object returned by \code{BOSS_modal()}, which must
#'   contain \code{design_points$x_original} and \code{objective_function}.
#' @param approach Character string specifying the method used for estimating the Hessian at the BO mode.
#' Must be one of \code{num.GP.refine}, \code{num.obj}, \code{num.GP}.
#' @param GP.refine.args List of additional arguments for the \code{num.GP.refine} method.
#' @param num.args List of additional arguments for \code{numDeriv::hessian()}.
#' @param tol Numeric tolerance for eigenvalues (default \code{1e-8}).
#' @param ... Further arguments passed to \code{numDeriv::hessian()}.
#'
#' @return The input \code{boss_result} with two new components:
#'   \itemize{
#'     \item \code{mode_hessian}: symmetric Hessian matrix at the mode.
#'     \item \code{mode}: the mode point in original input space.
#'   }
#'
#' @details
#' The function numerically computes the hessian at the mode obtained from a \code{boss} Object using three different strategies:
#'
#' 1. Default: Numerical hessian based on refined GP surrogate (\code{num.GP.refine}): given a specified accuracy level \code{GP.refine.args$eps},
#' the function first draws a circle with radius \code{GP.refine.args$eps/2} around the mode from the \code{boss} object and check how many \code{boss} design points
#' are within the circle. If the number is below \eqn{GP.refine.args$n_add}, additional uniformly distributed design points are added within the circle to reach \eqn{GP.refine.args$n_add} and
#' evaluate at these additional points with \code{boss_result$objective_function}. Combining all design points in the circle, estimate the Hessian at \code{boss} mode by
#' fitting another GP surrogate.
#' If \code{GP.refine.args$eps = NULL},
#' it is obtained by
#' \deqn{\epsilon = \frac{1}{4\sqrt{\text{Tr}(-H_{GP})}},}
#' where \eqn{H_{GP}} is the hessian at the \code{boss} mode computed from \code{numDeriv::hessian()} on the non-refined GP surrogate \code{boss_result$surrogate}.
#' If \code{GP.refine.args$n_add = NULL}, it is set to \code{(D+2)*(D+1)}.
#'
#' 2. Brute-force numerical hessian (\code{num.obj}): directly estimate the hessian at the \code{boss} mode via \code{numDeriv::hessian()} using the \code{boss_result$objective_function}.
#'
#' 3. Numerical hessian based on GP surrogate (\code{num.GP}): estimate the hessian at the \code{boss} mode via \code{numDeriv::hessian()} using the \code{boss_result$surrogate}.
#'
#' Note that \code{local.poly} balances between theoretical accuracy and computational budget. \code{num.obj} is the most computationally intense while \code{num.GP} is the cheapest, but does not have theoretical guarantee.
#' In addition, if there is noise in the evaluation of \code{boss_result$objective_function}, we recommend users to use \code{local.poly} or \code{num.GP}. Otherwise, \code{numDeriv::hessian()} based on \code{boss_result$objective_function} will be highly unstable.
#'
#' @importFrom numDeriv hessian
#' @importFrom lhs randomLHS
#' @export
update_hessian <- function(boss_result,
                           approach      = 'num.GP.refine',
                           GP.refine.args = list(eps = NULL, n_add = NULL),
                           num.args = list(method = 'Richardson', method.args = list()),
                           tol         = 1e-8,
                           ...) {
  if(!inherits(boss_result, 'boss_modal')){
    stop(paste0('No applicable method from update_hessian() to object of class ', class(boss_result)))
  }
  # extract original (unscaled) design matrix and responses
  xmat_orig <- boss_result$design_points$x_original
  yvec      <- boss_result$design_points$y
  D         <- ncol(xmat_orig)

  # locate current optimizer in original space
  mode_point <- xmat_orig[which.max(yvec), , drop = FALSE]

  if(!approach %in% c('num.GP.refine', 'num.obj', 'num.GP')){
    stop("Hessian estimation approach must be one of 'num.GP.refine', 'num.obj', 'num.GP'.")
  }
  else if(approach == 'num.obj'){
    # compute Hessian via numDeriv based on objective function
    H <- numDeriv::hessian(
      func        = boss_result$objective_function,
      x           = mode_point,
      method      = num.args$method,
      method.args = num.args$method.args,
      ...
    )
    # symmetrize
    H <- (H + t(H)) / 2
  }
  else if(approach == 'num.GP'){
    # compute Hessian via numDeriv based on GP surrogate
    if(!boss_result$quad) {
      xmat_s <- sweep(sweep(xmat_orig, 2, boss_result$lower), 2, boss_result$upper - boss_result$lower, FUN = '/')

      opt <- optim(runif(1, 0.01, 0.99), function(l)
        compute_like(length_scale = l,
                     y = yvec,
                     x = xmat_s,
                     quad = T,
                     D = D,
                     noise_var = boss_result$noise_var),
        control = list(maxit = 100), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
      length_scale <- opt$par

      signal_var <- compute_like(length_scale = NULL, y = yvec, x = xmat_s, quad = T, D = D)

      s_choice_cov <- cov_generator(length_scale = length_scale,
                                    signal_var = signal_var,
                                    nu = boss_result$gp_params$nu)

      s_surrogate <- function(xvalue, data_to_smooth) {
        xvalue_s <- (xvalue - boss_result$lower) / (boss_result$upper - boss_result$lower)
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue_s, ncol = D),
                   choice_cov = s_choice_cov, noise_var = boss_result$noise_var,
                   quad = T)$mean
      }

      s_design <- list(x = xmat_s, y = yvec)
      fn_s <- function(y) as.numeric(s_surrogate(xvalue = y, data_to_smooth = s_design))
    }
    else{
      fn_s <- boss_result$surrogate
    }

    H_crude <- numDeriv::hessian(
      func        = boss_result$surrogate,
      x           = mode_point,
      method      = num.args$method,
      method.args = num.args$method.args,
      ...
    )

    eps <- 0.25/sqrt((sum(diag(abs(H_crude)))))

    H <- optim(mode_point, fn_s, method = 'L-BFGS-B',
               control = list(maxit = 100, fnscale = -1), lower = mode_point - eps, upper = mode_point + eps,
               hessian = T)$hessian

    H <- (H + t(H))/2
  }
  else{
    # compute Hessian via locally weighted polynomial regression
    eps <- GP.refine.args$eps

    if(is.null(eps)){

      H_crude <- numDeriv::hessian(
        func        = boss_result$surrogate,
        x           = mode_point,
        method      = num.args$method,
        method.args = num.args$method.args,
        ...
      )

      eps <- 0.25/sqrt((sum(diag(abs(H_crude)))))
    }

    required_neighbors <- GP.refine.args$n_add
    if(is.null(required_neighbors)){
      required_neighbors <- (D+2)*(D+1)
    }

    dist <- sqrt(rowSums(sweep(xmat_orig, 2, mode_point)^2))
    nn_idx <- dist < eps/2
    xmat_neighbor <- xmat_orig[nn_idx, , drop = FALSE]
    xmat_not_neighbor <- xmat_orig[!nn_idx, , drop=FALSE]
    y_neighbor <- yvec[nn_idx]
    y_not_neighbor <- yvec[!nn_idx]

    if(sum(nn_idx) < required_neighbors){

      U_dir <- matrix(runif(required_neighbors*D), ncol = D)
      Z <- qnorm(U_dir)                          # Map to N(0,1)
      Z <- Z / sqrt(rowSums(Z^2))               # Normalize to unit sphere

      # LHS for radius component
      r <- lhs::randomLHS(required_neighbors, 1)[,1]^(1 / D)  # Ensure shape is numeric vector
      X_add <- Z * (eps / 2 * r)     # Scale by radius

      # Shift to center at mode_point
      X_add <- sweep(X_add, 2, mode_point, FUN = "+")

      dists <- apply(X_add, 1, function(xa) {
        min(sqrt(rowSums((t(t(xmat_orig) - xa))^2)))
      })

      # Get indices of the n closest points to any xmat_original point
      remove_idx <- order(dists)[1:sum(nn_idx)]

      # Remove them
      X_add[-remove_idx, , drop = FALSE]

      # Evaluate objective function at new points
      y_add <- apply(X_add, 1, boss_result$objective_function)

      # Combine neighbors and new LHS-sampled points
      X_local <- rbind(xmat_neighbor, X_add)
      y_local <- c(y_neighbor, y_add)
    }
    else{
      X_local <- xmat_neighbor
      y_local <- y_neighbor
    }

    X_refined <- rbind(X_local, xmat_not_neighbor)
    y_refined <- c(y_local, y_not_neighbor)

    X_refined_s <- sweep(sweep(X_refined, 2, boss_result$lower), 2, boss_result$upper - boss_result$lower, FUN = '/')

    opt <- optim(runif(1, 0.01, 0.99), function(l)
      compute_like(length_scale = l,
                   y = y_refined,
                   x = X_refined_s,
                   quad = T,
                   D = D,
                   noise_var = boss_result$noise_var),
      control = list(maxit = 100), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
    length_scale <- opt$par

    signal_var <- compute_like(length_scale = NULL, y = y_refined, x = X_refined_s,
                               quad = T, D = D)

    refined_choice_cov <- cov_generator(length_scale = length_scale,
                                        signal_var = signal_var,
                                        nu = boss_result$gp_params$nu)

    refined_surrogate <- function(xvalue, data_to_smooth) {
      xvalue_s <- (xvalue - boss_result$lower) / (boss_result$upper - boss_result$lower)
      predict_gp(data = data_to_smooth, x_pred = matrix(xvalue_s, ncol = D),
                 choice_cov = refined_choice_cov, noise_var = boss_result$noise_var,
                 quad = T)$mean
    }

    refined_design <- list(x = X_refined_s, y = y_refined)
    fn_refined <- function(y) as.numeric(refined_surrogate(xvalue = y, data_to_smooth = refined_design))

    H <- optim(mode_point, fn_refined, method = 'L-BFGS-B',
               control = list(maxit = 100, fnscale = -1), lower = mode_point - eps, upper = mode_point + eps,
               hessian = T)$hessian

    H <- (H + t(H))/2
  }

  # check negative semi-definiteness: all eigenvalues <= tol
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev > tol)) {
    stop(sprintf(
      "Computed Hessian is not negative semi-definite (max eigenvalue = %g > %g)",
      max(ev), tol
    ))
  }

  boss_result$mode_hessian <- H
  boss_result$mode         <- mode_point[ , , drop = TRUE]
  return(boss_result)
}




#' Compute Essential Support for a \code{boss_moal} Object
#'
#' Compute and store the 1-\eqn{\alpha} level set (ellipsoidal support) defined by the
#' Hessian at the mode in a \code{boss_modal} object.
#'
#' @param boss_result A S3 \code{boss_modal} object (output of \code{BOSS_modal}), which must contain:
#'   \describe{
#'     \item{\code{design_points$x}}{Scaled design matrix.}
#'     \item{\code{design_points$y}}{Responses.}
#'     \item{\code{mode_hessian}}{Hessian matrix at the mode (from \code{update_hessian}).}
#'   }
#' @param alpha Numeric in [0,1], significance level (default 0.05) for the chi-squared cutoff.
#'
#' @return The input \code{boss_result} with a new or overwritten component
#'   \code{essential_support}, a list containing:
#'   \describe{
#'     \item{\code{center}}{Mode (vector of length D).}
#'     \item{\code{H}}{Hessian matrix at the mode.}
#'     \item{\code{chi2_radius}}{Scalar \eqn{c = \chi^2_{D,1-\alpha}}.}
#'     \item{\code{eig_vectors}}{Eigenvectors of \code{H}.}
#'     \item{\code{eig_values}}{Eigenvalues of \code{H}.}
#'   }
#'
#' @examples
#' # Assuming `br` is a boss object with mode_hessian (Not run)
#' # br <- compute_essential_support(br, alpha = 0.05)
#' # str(br$essential_support)
#'
#' @importFrom stats qchisq
#' @export
compute_essential_support <- function(boss_result, alpha = 0.05) {
  if(!inherits(boss_result, 'boss_modal')){
    stop(paste0('No applicable method from compute_essential_support() to object of class ', class(boss_result)))
  }
  if (!is.null(boss_result$essential_support)) {
    warning("Overwriting existing essential_support in boss_result")
  }
  if (is.null(boss_result$mode_hessian)) {
    stop("mode_hessian not found. Please run update_hessian() first.")
  }

  H <- boss_result$mode_hessian
  # Determine the mode (max response)
  center <- boss_result$mode
  D      <- length(center)

  # Chi-square cutoff
  chi2_radius <- stats::qchisq(1 - alpha, df = D)

  # Eigen-decomposition
  eig <- eigen(H)

  # Store essential support
  boss_result$essential_support <- list(
    center      = center,
    H           = H,
    chi2_radius = chi2_radius,
    eig_vectors = eig$vectors,
    eig_values  = eig$values
  )

  return(boss_result)
}

