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
#' @importFrom fields rdist
#' @importFrom randtoolbox sobol
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
#' \deqn{\sup_{x\in\Omega}\min_j \|x - x_j\|}
#' by:
#' \enumerate{
#' \item Uniformly sample inside the ellipsoid \eqn{\Omega},
#' \item Truncating all candidates to lie within the original \code{lower}/\code{upper} box,
#' \item Computing for each candidate its nearest-design-point distance from the uniform sample,
#' \item Taking the maximum of these minima.
#' }
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support}, \code{essential_design_points$x_original},
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

#' Fill-in essential support with a Sobol-sequence
#'
#' Given a \code{boss} object with \code{essential_support} and \code{essential_design_points$x_original}, this function will
#' approximately fill-in with a Sobol-sequence in the essential support and filter against existing design points:
#' \enumerate{
#' \item Based on the required fill-in distance \code{h}, estimate the required number of quasi-uniform design points by
#'   \deqn{
#'     n \;=\; \max\Bigl\{\texttt{n\_sample\_max},\;
#'                      2\,\mathrm{Vol}(\text{Original box})\;\big/\;h^{D}\Bigr\}.
#'   }
#' \item Generate a Sobol-sequence within the original box.
#' \item Truncate the Sobol-sequence outside the ellipsoid.
#' \item Compute the distance from existing design points to the Sobol-sequence and remove Sobol-candidates that are within \code{h} of existing points.
#' \item Add in a maximum of \code{max_add} number selected fill-in Sobol-candidates with optimal covering over the ellipsoid ensured.
#' \item Recompute new fill-in based on the selected Sobol and existing design points.
#' }
#'
#' If the candidate pool is over \code{max_add}, or if \code{fill_in} cannot be reduced below \code{h},
#' a warning is issued.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support} and \code{essential_design_points$x_original},
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


#' Update a \code{boss} Object After Adding New Essential Fill-in Design Points
#'
#' Given a \code{boss} object that already has an \code{essential_design_points} list,
#' this function does the following:
#' \enumerate{
#'   \item Fills in any missing \code{NA} responses by evaluating the objective function.
#'   \item Recomputes the mode and its Hessian based on the updated design.
#'   \item Re-fits the GP hyperparameters (length-scale and signal variance) using MLE.
#'   \item Constructs and stores an efficient \code{surrogate} closure for fast prediction,
#'         using precomputed posterior weights (and GLS mean model if \code{quad=TRUE}).
#' }
#'
#'
#' @param boss_result A \code{boss} object with fields:
#'   \describe{
#'     \item{\code{essential_design_points}}{List with \code{x_original}, \code{y}, \code{x}.}
#'     \item{\code{lower}, \code{upper}}{Bounds for scaling inputs.}
#'     \item{\code{objective_function}}{Callable to evaluate new design points.}
#'     \item{\code{gp_params}}{Includes \code{noise_var}, \code{nu}, \code{quad} flag, etc.}
#'   }
#'
#' @return The updated \code{boss_result}, with:
#'   \describe{
#'     \item{\code{essential_design_points$y}}{Updated with evaluated values.}
#'     \item{\code{mode}, \code{mode_hessian}}{Recomputed for the surrogate.}
#'     \item{\code{gp_params}}{Updated hyperparameters.}
#'     \item{\code{surrogate}}{A closure for predicting at new inputs.}
#'   }
#'
#' @export
update_boss <- function(boss_result) {
  ## 1) Fill missing y's in essential_design_points
  ed <- boss_result$essential_design_points
  nas <- which(is.na(ed$y))
  if (length(nas) > 0) {
    for (i in nas) {
      ed$y[i] <- boss_result$objective_function(ed$x_original[i, ])
    }
    boss_result$essential_design_points$y <- ed$y
  }

  ## 2) Temporarily override design_points
  old_dp <- boss_result$design_points
  boss_result$design_points <- list(
    x_original = ed$x_original,
    x = sweep(ed$x_original, 2, boss_result$lower, "-") /
      (boss_result$upper - boss_result$lower),
    y = ed$y
  )

  ## 3) Update mode and Hessian
  boss_result <- update_hessian(boss_result)

  ## 4) Refit GP hyperparameters
  x_scaled <- boss_result$design_points$x
  y_obs    <- boss_result$design_points$y - mean(boss_result$design_points$y)

  opt <- optim(
    runif(1, 0.01, 0.99),
    function(l) compute_like(
      length_scale = l,
      y = y_obs,
      x = x_scaled,
      D = boss_result$D,
      nu = boss_result$gp_params$nu,
      quad = boss_result$quad,
      noise_var = boss_result$noise_var
    ),
    control = list(maxit = 1000),
    lower = 0.01,
    upper = 0.99,
    method = 'L-BFGS-B'
  )

  signal_var <- compute_like(
    length_scale = NULL,
    y = y_obs,
    x = x_scaled,
    D = boss_result$D,
    quad = boss_result$quad,
  )

  boss_result$gp_params$length_scale <- opt$par
  boss_result$gp_params$signal_var   <- signal_var

  ## 5) Precompute for prediction
  intern <- predict_gp_internal(
    data = list(
      x = boss_result$essential_design_points$x,
      y = (boss_result$essential_design_points$y - mean(boss_result$essential_design_points$y))
    ),
    noise_var = boss_result$noise_var,
    choice_cov = cov_generator(
      length_scale = boss_result$gp_params$length_scale,
      nu           = boss_result$gp_params$nu,
      signal_var   = boss_result$gp_params$signal_var
    ),
    quad = boss_result$quad
  )

  ## 6) Build surrogate closure
  boss_result$surrogate <- function(xnew){
    # Scale
    xnew <- matrix(xnew, ncol = boss_result$D)
    x_s <- sweep(xnew, 2, boss_result$lower, "-") / (boss_result$upper - boss_result$lower)
    obtain_mean_internal(
      xnew = x_s, intern = intern,
      covfn = cov_generator(
        length_scale = boss_result$gp_params$length_scale,
        nu           = boss_result$gp_params$nu,
        signal_var   = boss_result$gp_params$signal_var),
      D = boss_result$D) + mean(boss_result$essential_design_points$y)
  }

  ## 7) Restore original design_points
  boss_result$design_points <- old_dp

  return(boss_result)
}









