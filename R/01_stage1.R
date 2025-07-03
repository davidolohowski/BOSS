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
#' @param nu GP covariance smoothness. \eqn{\nu = \infty} represents squared-exponential kernel; \eqn{0 < \nu < \infty} represents the Matern kernel with smoothness \eqn{\lceil\nu\rceil - 1}.
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
#' @return A "boss" object (S3 class \code{boss}), which is a list containing:
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
  signal_var <- var(yvec)

  # Initial optimization for lengthscale.
  if (!is.null(opt.lengthscale.grid)) {
    length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
    like_vec <- sapply(length_scale_vec, function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var))
    max_idx <- which.max(like_vec)
    length_scale <- length_scale_vec[max_idx]
    lik <- like_vec[max_idx]
  } else {
    opt <- optim(runif(1, 0.01, 0.99), function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var),
      control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
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
      signal_var <- var(newdata$y)
      if(!is.null(opt.lengthscale.grid)) {
        length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
        like_vec <- sapply(length_scale_vec, function(l)
          -compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                        signal_var = signal_var, noise_var = noise_var))
        max_idx <- which.max(like_vec)
        length_scale <- length_scale_vec[max_idx]
        lik <- like_vec[max_idx]
      } else {
        opt <- optim(runif(1, 0.01, 0.99), function(l)
          compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                       signal_var = signal_var, noise_var = noise_var),
          control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
        length_scale <- opt$par
        lik <- opt$value
      }
      if(verbose == 3) {
        print(paste("The new length.scale:", length_scale))
        print(paste("The new signal_var:", signal_var))
      }
      choice_cov <- cov_generator(length_scale = length_scale, signal_var = signal_var, nu = nu)
    }

    # Optimize the acquisition function (UCB).
    if(verbose == 3) print("Maximize Acquisition Function")
    if(is.null(opt.grid)){
      # Multi-start local optimization.
      initialize_UCB <- matrix(runif(D*optim.n, rep(0, D), rep((1), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            quad = quad,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
        control = list(maxit = 100),
        lower = rep(0, D), upper = rep((1), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
             quad = quad,
             cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
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

    yvec <- c(yvec + rel, (y_new))
    rel <- mean(yvec)
    yvec <- yvec - rel

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

  surrogate_fn <- function(xvalue) {
    xvalue_transform <- (xvalue - lower) / (upper - lower)
    predict_gp(
      data    = list(x = xmat_trans, y = (yvec)),
      x_pred  = matrix(xvalue_transform, ncol = D),
      choice_cov = choice_cov,
      noise_var  = noise_var,
      quad       = quad
    )$mean + rel
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
    class = "boss")

  return(boss_result)
}




#' Update Mode Hessian for a \code{boss} Object (with PSD check)
#'
#' Compute the numerical Hessian of the objective at the current best point
#' (in original input space), verify it is negative semi-definite (NSD), and
#' store it back into the \code{boss_result}.
#'
#' @param boss_result A S3 \code{boss} object returned by \code{BOSS_modal()}, which must
#'   contain \code{design_points$x_original} and \code{objective_function}.
#' @param method Character string specifying the finite-difference method
#'   passed to \code{numDeriv::hessian} (default \code{"Richardson"}).
#' @param method.args List of additional arguments for the chosen method.
#' @param tol Numeric tolerance for eigenvalues (default \code{1e-8}).
#' @param ... Further arguments passed to \code{numDeriv::hessian}.
#'
#' @return The input \code{boss_result} with two new components:
#'   \itemize{
#'     \item \code{mode_hessian}: symmetric Hessian matrix at the mode.
#'     \item \code{mode}: the mode point in original input space.
#'   }
#' @importFrom numDeriv hessian
#' @export
update_hessian <- function(boss_result,
                           method      = "Richardson",
                           method.args = list(),
                           tol         = 1e-8,
                           ...) {
  # extract original (unscaled) design matrix and responses
  xmat_orig <- boss_result$design_points$x_original
  yvec      <- boss_result$design_points$y
  D         <- ncol(xmat_orig)

  # locate current optimizer in original space
  mode_point <- xmat_orig[which.max(yvec), , drop = TRUE]

  # compute Hessian via numDeriv
  H <- numDeriv::hessian(
    func        = boss_result$objective_function,
    x           = mode_point,
    method      = method,
    method.args = method.args,
    ...
  )
  # symmetrize
  H <- (H + t(H)) / 2

  # check negative semi-definiteness: all eigenvalues <= tol
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev > tol)) {
    stop(sprintf(
      "Computed Hessian is not negative semi-definite (max eigenvalue = %g > %g)",
      max(ev), tol
    ))
  }

  boss_result$mode_hessian <- H
  boss_result$mode         <- mode_point
  return(boss_result)
}




#' Compute Essential Support for a \code{boss} Object
#'
#' Compute and store the 1-\eqn{\alpha} level set (ellipsoidal support) defined by the
#' Hessian at the mode in a \code{boss} object.
#'
#' @param boss_result A S3 \code{boss} object (output of \code{BOSS_modal}), which must contain:
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

