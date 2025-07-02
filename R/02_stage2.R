#' Extract Design Points within the essential support for a \code{boss} Object
#'
#' Given a \code{boss} object with \code{essential_support} computed (in original input space),
#' builds \code{essential_design_points} that contains all design points inside the ellipsoid.
#'
#' @param boss_result A \code{boss} object containing non-NULL
#'   \code{essential_support}, and with \code{lower} and \code{upper} fields.
#'
#' @return The same \code{boss_result}, with
#'   \code{essential_design_points} list:
#'   \describe{
#'     \item{\code{x}}{Scaled points inside support (matrix).}
#'     \item{\code{x_original}}{Original-space points inside support (matrix).}
#'     \item{\code{y}}{Responses: real values for design points inside support.}
#'   }
#'
#' @export
construct_essential_designs <- function(boss_result) {
  # prerequisites
  ess <- boss_result$essential_support
  if (is.null(ess)) stop("essential_support must be computed first")
  if (is.null(boss_result$lower) || is.null(boss_result$upper)) {
    stop("boss_result must contain `lower` and `upper` for truncation")
  }

  # support parameters (original space)
  mu    <- ess$center
  H     <- -ess$H            # flip to PSD
  cval  <- ess$chi2_radius
  D     <- length(mu)

  # original and scaled design
  Xs <- boss_result$design_points$x
  Xo <- boss_result$design_points$x_original
  Y  <- boss_result$design_points$y

  # interior points
  diffs <- t(t(Xo) - mu)
  quad  <- rowSums((diffs %*% H) * diffs)
  idx   <- which(quad <= cval)
  Xs_int  <- Xs[idx,, drop=FALSE]
  Xo_int  <- Xo[idx,, drop=FALSE]
  Y_int   <- Y[idx]

  boss_result$essential_design_points <- list(
    x          = Xs_int,
    x_original = Xo_int,
    y          = Y_int
  )
  return(boss_result)
}


#' Approximate Fill-In Distance via Uniform Sampling (with Truncation)
#'
#' Estimate the fill-in distance
#' \(\sup_{x\in\Omega}\min_j \|x - x_j\|\)
#' by:
#' 1. Uniformly sample inside the ellipsoid Ω,
#' 2. Truncating all candidates to lie within the original \code{lower}/\code{upper} box,
#' 3. Computing for each candidate its nearest-design-point distance from the uniform sample,
#' 4. Taking the maximum of these minima.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support}, \code{essential_design_points\$x_original},
#'   and containing \code{lower} and \code{upper}.
#' @param n_samples Integer; number of unifrom samples in the essential support to estimate the fill-in (default 10000).
#'
#' @return The input \code{boss_result}, with \code{fill_in} updated.
#'
#' @importFrom stats dist
#' @export
compute_fill_in <- function(boss_result, n_samples = 10000) {
  ess <- boss_result$essential_support
  ed <- boss_result$essential_design_points
  if (is.null(ess) || is.null(ed$x_original)) {
    stop("Please run compute_essential_support() and construct_essential_designs() first.")
  }

  center <- ess$center
  Hinv <- -ess$H
  cval <- ess$chi2_radius
  Xd <- as.matrix(ed$x_original)
  if (is.vector(Xd)) {
    Xd <- matrix(Xd, ncol = 1)
  }
  lower <- boss_result$lower
  upper <- boss_result$upper
  D <- length(center)

  Ainv <- solve(Hinv)

  sample_ellipsoid <- function(n, center, Ainv, cval) {
    unit_pts <- matrix(runif(n * D), ncol = D)
    z <- qnorm(unit_pts)
    z <- z / sqrt(rowSums(z^2))
    r <- runif(n)^(1 / D)
    z <- z * r

    L <- chol(Ainv)
    scaled <- z %*% L * sqrt(cval)
    sweep(scaled, 2, center, "+")
  }

  test_pts <- sample_ellipsoid(n_samples, center, Ainv, cval)
  test_pts <- pmin(pmax(test_pts, matrix(lower, nrow = n_samples, ncol = D, byrow = TRUE)),
                   matrix(upper, nrow = n_samples, ncol = D, byrow = TRUE))

  dmat <- fields::rdist(test_pts, Xd)
  min_d <- apply(dmat, 1, min)

  boss_result$fill_in <- max(min_d, na.rm = TRUE)
  return(boss_result)
}

#' Approximately fill-in by supplying a Sobol-sequence in the essential support and filter against existing design points.
#'
#' 1. Based on the required fill-in distance \code{h}, compute the required number of quasi-uniform design points by
#'    \(n = \max\{n_sample_max, 2*\text{Vol}(\text{Original box})/h^D\}\).
#' 2. Generate a Sobol-sequence within the original box.
#' 3. Truncate the Sobol-sequence within the ellipsoid.
#' 4. Compute the distance from existing design points to the Sobol-sequence and remove Sobol-candidates that are within h of existing points.
#' 5. Add in a maximum of \code{max_add} number selected fill-in Sobol-candidates to ensure optimal covering over the ellipsoid.
#' 6. Recompute new fill-in based on the selected Sobol and existing design points.
#'
#' If the candidate pool is over \code{max_add}, or if \code{fill_in} cannot be reduced below \code{h},
#' a warning is issued.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support} and \code{essential_design_points\$x_original},
#'   and containing \code{lower} and \code{upper}.
#' @param h Numeric; target fill-in distance (>= 0).
#' @param max_add Integer; maximum number of points to add (default 100).
#' @param n_sample_max Integer; maximum number of Sobol-sequence candidates to be added (default 10000).
#' @param verbose Integer in \{0,1\}; 0=silent, 1=summary.
#'
#' @return The updated \code{boss_result}, with its
#'   \code{essential_design_points} and \code{fill_in} refreshed.
#'
#' @importFrom stats dist
#' @export
fill_in <- function(boss_result, h, max_add = 100, n_sample_max = 10000, verbose = 0) {
  ess <- boss_result$essential_support
  ed <- boss_result$essential_design_points
  if (is.null(ess) || is.null(ed$x_original)){
    stop('Please run compute_essential_support() and construct_essential_designs() first')
  }
  if (!is.numeric(boss_result$fill_in)){
    stop('Please run compute_fill_in() before calling fill_in()')
  }
  if (!is.numeric(h) || length(h) != 1 || h <= 0) {
    stop("`h` must be a positive numeric scalar")
  }
  if(!verbose %in% 0:1){
    stop("`verbose` must be an integer in {0,1}")
  }

  center <- ess$center
  Hinv <- -ess$H
  cval <- ess$chi2_radius
  lower <- boss_result$lower
  upper <- boss_result$upper
  D <- length(center)

  Xo <- as.matrix(ed$x_original)
  if (D == 1) {
    Xo <- matrix(Xo, ncol = 1)
  }
  Xs <- as.matrix(ed$x)
  if (D == 1) {
    Xs <- matrix(Xs, ncol = 1)
  }
  Y <- ed$y

  # Estimate number of samples
  box_vol <- prod(upper - lower)
  n_samples <- min(ceiling(2*box_vol / h^D), n_sample_max)
  if (ceiling(2*box_vol / h^D) >= n_sample_max){
    warning('Required number of fill-in points to reach h is greater than n_sample_max. Consider increase n_sample_max.')
  }

  # Sobol sampling in the box
  sobol_unit <- as.matrix(randtoolbox::sobol(n = n_samples, dim = D, scrambling = 0))
  sobol_unit <- matrix(sobol_unit, ncol = D)
  sobol_box <- sweep(sobol_unit, 2, upper - lower, "*") + matrix(lower, nrow = n_samples, ncol = D, byrow = TRUE)

  # Filter by ellipsoid membership
  diffs_cand <- sweep(sobol_box, 2, center)
  quad_cand <- rowSums((diffs_cand %*% Hinv) * diffs_cand)
  sobol_inside <- sobol_box[quad_cand <= cval, , drop = FALSE]

  # Distance to existing points
  dmat <- if (nrow(Xo) > 0) fields::rdist(sobol_inside, Xo) else matrix(Inf, nrow = nrow(sobol_inside), ncol = 1)
  min_dist <- apply(dmat, 1, min)

  # Filter by fill-in distance
  keep_idx <- which(min_dist > h)
  candidates <- sobol_inside[keep_idx, , drop = FALSE]

  if (nrow(candidates) == 0) {
    if (verbose >= 1) {
      cat("No additional points needed — fill-in already satisfied.\n")
    }
    boss_result$essential_design_points$x_original <- Xo
    boss_result$essential_design_points$x <- Xs
    boss_result$essential_design_points$y <- Y
    return(boss_result)
  }

  n_add <- min(nrow(candidates), max_add)
  if(nrow(candidates) >= max_add){
    warning('Number of points to be added is greater than max_add. Required h may not be achieved.')
  }
  new_pts <- candidates[seq_len(n_add), , drop = FALSE]
  new_scaled <- sweep(new_pts, 2, lower, "-") / (upper - lower)

  Xo <- rbind(Xo, new_pts)
  Xs <- rbind(Xs, new_scaled)
  Y <- c(Y, rep(NA_real_, n_add))

  boss_result$essential_design_points$x_original <- Xo
  boss_result$essential_design_points$x <- Xs
  boss_result$essential_design_points$y <- Y

  boss_result <- compute_fill_in(boss_result, n_samples = 100*nrow(Xo))

  if (verbose >= 1) {
    cat("fill in: added", n_add, "point(s).\n")
  }
  if (boss_result$fill_in > h){
    warning(sprintf('Updated fill-in distance is (%.6f) > target h (%.6f); Adjust your expectation by either increasing max_add and n_sample_max or increasing h.',
                    boss_result$fill_in, h))
  }

  return(boss_result)
}




#' Update a \code{boss} Object After Filling Essential Design Points
#'
#' Given a \code{boss_result} with \code{essential_design_points} already constructed,
#' this function will:
#' 1. Fill in any missing responses (\code{NA}) by evaluating the objective.
#' 2. Recompute the mode and Hessian at the mode.
#' 3. Re‐fit the GP hyperparameters (length‐scale and signal variance).
#' 4. Rebuild the surrogate function closure with the updated design.
#'
#' @param boss_result A \code{boss} object with
#'   \code{essential_design_points} (list with \code{x_original}, \code{y}),
#'   and containing \code{lower}, \code{upper}, and \code{objective_function}.
#' @return The updated \code{boss_result}, with fields
#'   \code{essential_design_points\$y}, \code{mode}, \code{mode_hessian},
#'   \code{gp_params}, and \code{surrogate} refreshed.
#' @export
update_boss <- function(boss_result) {
  ## 1) Fill missing y's
  ed <- boss_result$essential_design_points
  nas <- which(is.na(ed$y))
  if (length(nas) > 0) {
    for (i in nas) {
      ed$y[i] <- boss_result$objective_function(ed$x_original[i, ])
    }
    boss_result$essential_design_points$y <- ed$y
  }

  ## 2) Update mode & Hessian based on essential_design_points
  # Temporarily override design_points to the essential ones
  old_dp_xo <- boss_result$design_points$x_original
  old_dp_x  <- boss_result$design_points$x
  old_dp_y  <- boss_result$design_points$y

  boss_result$design_points$x_original <- ed$x_original
  # recompute scaled x
  boss_result$design_points$x <- sweep(ed$x_original, 2,
                                       boss_result$lower, "-") /
    (boss_result$upper - boss_result$lower)
  boss_result$design_points$y <- ed$y

  # update mode and Hessian
  boss_result <- update_hessian(boss_result)

  ## 3) Re-fit GP hyperparameters on the essential design
  # use compute_like + optim
  x_scaled <- boss_result$design_points$x
  y_centered <- boss_result$design_points$y - mean(boss_result$design_points$y)
  signal_var <- var(y_centered)

  opt <- optim(runif(1, 0.01, 0.99), function(l)
    compute_like(length_scale = l, y = y_centered, x = x_scaled, D = boss_result$D,
                 nu = boss_result$gp_params$nu, quad = boss_result$quad,
                 signal_var = signal_var, noise_var = boss_result$noise_var),
    control = list(maxit = 1000), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')

  length_scale <- opt$par
  boss_result$gp_params$length_scale <- length_scale
  boss_result$gp_params$signal_var   <- signal_var

  ## 4) Rebuild surrogate closure
  boss_result$surrogate <- function(xnew) {
    xnew_transformed <- (xnew - boss_result$lower)/
      (boss_result$upper - boss_result$lower)
    predict_gp(data        = list(x = boss_result$essential_design_points$x,
                                  y = boss_result$essential_design_points$y - mean(boss_result$essential_design_points$y)),
               x_pred      = matrix(xnew_transformed, ncol = boss_result$D),
               choice_cov  = cov_generator(
                 length_scale = boss_result$gp_params$length_scale,
                 nu = boss_result$gp_params$nu,
                 signal_var   = boss_result$gp_params$signal_var),
               noise_var   = boss_result$noise_var,
               quad        = boss_result$quad)$mean + mean(boss_result$essential_design_points$y)
  }

  ## restore original design_points
  boss_result$design_points$x_original <- old_dp_xo
  boss_result$design_points$x          <- old_dp_x
  boss_result$design_points$y          <- old_dp_y

  return(boss_result)
}









