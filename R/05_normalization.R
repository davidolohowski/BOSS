#' @title Add normalized posterior to boss object
#' @description
#' Given a boss object with a surrogate log-density function and precomputed essential support,
#' this function estimates the normalization constant (either numerically or via Monte Carlo)
#' and adds a normalized posterior density function as a new attribute in the boss object.
#'
#' It assumes that boss$essential_support$H is the Hessian of the log-posterior (negative definite at the mode).
#' All internal computations automatically negate it to get the positive-definite precision approximation.
#'
#' @param boss A boss object with at least $surrogate (log-density function), $D (dimension),
#'   and $essential_support (computed by compute_essential_support()).
#' @param int_method Integration method. One of "numeric" (1D only) or "mc".
#' @param nsamples Number of samples to use for Monte Carlo integration (default 10000).
#' @param ... Additional arguments passed to integrate() in numeric mode.
#'
#' @return The updated boss object with $normalized_posterior attribute.
#' @export
get_normalized_posterior <- function(
    boss,
    int_method = c("numeric", "mc"),
    nsamples = 10000,
    ...
) {
  int_method <- match.arg(int_method)

  if (is.null(boss$essential_support)) {
    stop("boss$essential_support is missing. Please run compute_essential_support() first.")
  }

  if (is.null(boss$D)) {
    stop("boss$D (dimension) is missing in the boss object.")
  }

  log_density <- boss$surrogate

  if (int_method == "numeric") {
    if (boss$D > 1) {
      stop("Numeric integration is only supported for 1D problems. Please use int_method = 'mc' for D > 1.")
    }

    bounds <- .extract_numeric_bounds_from_ellipsoid(boss$essential_support)
    lower <- bounds[1]
    upper <- bounds[2]

    unnormalized_density <- function(x) exp(log_density(x))
    normalization_constant <- integrate(unnormalized_density, lower, upper, ...)$value

    boss$normalized_posterior <- function(x) {
      exp(log_density(x)) / normalization_constant
    }

  } else if (int_method == "mc") {
    samples <- .ellipsoid_sampler(nsamples, boss$essential_support)
    log_unnormalized <- apply(samples, 1, log_density)
    unnormalized <- exp(log_unnormalized)
    avg <- mean(unnormalized)

    vol <- .ellipsoid_volume(
      boss$essential_support$H,
      boss$essential_support$chi2_radius
    )

    Z_hat <- avg * vol

    boss$normalized_posterior <- function(x) {
      exp(log_density(x)) / Z_hat
    }
  }

  return(boss)
}

#' @keywords internal
.extract_numeric_bounds_from_ellipsoid <- function(essential_support) {
  if (length(essential_support$center) != 1) {
    stop("Numeric integration currently only supported for 1D problems.")
  }

  center <- essential_support$center
  H_posdef <- - essential_support$H
  chi2_radius <- essential_support$chi2_radius

  if (H_posdef <= 0) {
    stop("Hessian must be negative definite at the mode.")
  }

  radius <- sqrt(chi2_radius / H_posdef)
  c(center - radius, center + radius)
}


#' @keywords internal
.ellipsoid_sampler <- function(n, essential_support) {
  D <- length(essential_support$center)
  center <- essential_support$center
  H_posdef <- - essential_support$H
  eig <- eigen(H_posdef, symmetric=TRUE)
  eigvecs <- eig$vectors
  eigvals <- eig$values
  chi2_radius <- essential_support$chi2_radius

  if (any(eigvals <= 0)) {
    stop("Negative or zero eigenvalues in (-H). Hessian must be negative definite at the mode.")
  }

  # robust scaling
  if (D == 1) {
    scaling <- matrix(1 / sqrt(eigvals), nrow = 1, ncol = 1)
  } else {
    scaling <- eigvecs %*% diag(1 / sqrt(eigvals))
  }

  samples <- matrix(NA_real_, nrow = n, ncol = D)
  for (i in seq_len(n)) {
    u <- rnorm(D)
    u <- u / sqrt(sum(u^2))
    r <- runif(1)^(1/D)
    u <- u * r * sqrt(chi2_radius)
    samples[i, ] <- as.numeric(center + scaling %*% u)
  }

  return(samples)
}



#' @keywords internal
.ellipsoid_volume <- function(H_negdef, chi2_radius) {
  H_posdef <- - H_negdef
  D <- nrow(H_posdef)
  eigvals <- eigen(H_posdef, symmetric=TRUE, only.values=TRUE)$values
  if (any(eigvals <= 0)) {
    stop("Negative or zero eigenvalues in (-H). Hessian must be negative definite at the mode.")
  }

  volume_unit_ball <- pi^(D / 2) / gamma(D / 2 + 1)
  volume_unit_ball * (chi2_radius)^(D / 2) / sqrt(det(H_posdef))
}
