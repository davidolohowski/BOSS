#' @title Add normalized posterior to boss object
#' @description
#' Given a boss object with a surrogate log-density function and precomputed essential support,
#' this function estimates the normalization constant (via numeric, Monte Carlo, or AGHQ integration)
#' and adds a normalized posterior density function as a new attribute in the boss object.
#'
#' It assumes that boss$essential_support$H is the Hessian of the log-posterior (negative definite at the mode).
#' All internal computations automatically negate it to get the positive-definite precision approximation.
#'
#' @param boss A boss object with at least $surrogate (log-density function), $D (dimension),
#'   and $essential_support (computed by compute_essential_support()).
#' @param int_method Integration method. One of "numeric" (1D only), "mc", or "aghq".
#' @param nsamples Number of samples to use for Monte Carlo integration (default 10000).
#' @param aghq_K Number of quadrature points for AGHQ (default 4).
#' @param aghq_startingvalue Starting value in transformed y-space for AGHQ (optional).
#' @param truncate_posterior Whether to truncate the posterior to have zero mass outside the essential support.
#' @param ... Additional arguments passed to `integrate()` (for "numeric") or `aghq_bounded()` (for "aghq").
#'
#' @return The updated boss object with $normalized_posterior attribute.
#' @export
get_normalized_posterior <- function(
    boss,
    int_method = c("numeric", "mc", "aghq"),
    nsamples = 10000,
    aghq_K = 4,
    aghq_startingvalue = NULL,
    truncate_posterior = TRUE,
    ...
) {
  int_method <- match.arg(int_method)

  if (is.null(boss$essential_support)) {
    stop("boss$essential_support is missing. Please run compute_essential_support() first.")
  }

  if (is.null(boss$D)) {
    stop("boss$D (dimension) is missing in the boss object.")
  }

  D <- boss$D

  raw_log_density <- boss$surrogate
  center <- boss$essential_support$center
  H_negdef <- boss$essential_support$H
  chi2_radius <- boss$essential_support$chi2_radius
  log_density <- raw_log_density
  if (truncate_posterior){
    log_density <- function(x) {
      if (is.null(dim(x))) {
        if (length(x) == D) {
          # Single observation
          delta <- x - center
          qf <- -sum(delta * (H_negdef %*% delta))
          if (qf > chi2_radius) return(-Inf)
          return(raw_log_density(x))
        } else if (D == 1) {
          # Batch of 1D inputs
          out <- numeric(length(x))
          for (i in seq_along(x)) {
            xi <- x[i]
            delta <- xi - center
            qf <- -sum(delta * (H_negdef %*% delta))
            out[i] <- if (qf > chi2_radius) -Inf else raw_log_density(xi)
          }
          return(out)
        } else {
          stop("Input x must be of length D or matrix with D columns.")
        }
      } else {
        # x is a matrix: each row is a point
        out <- numeric(nrow(x))
        for (i in seq_len(nrow(x))) {
          xi <- x[i, ]
          delta <- xi - center
          qf <- -sum(delta * (H_negdef %*% delta))
          out[i] <- if (qf > chi2_radius) -Inf else raw_log_density(xi)
        }
        return(out)
      }
    }
  }


  if (int_method == "numeric") {
    if (D > 1) stop("Numeric integration is only supported for 1D problems.")
    bounds <- .extract_numeric_bounds_from_ellipsoid(boss$essential_support)
    lower <- bounds[1]
    upper <- bounds[2]

    normalization_constant <- integrate(function(x) exp(log_density(x)), lower, upper, ...)$value

    boss$normalized_posterior <- function(x) {
      exp(log_density(x)) / normalization_constant
    }

  } else if (int_method == "mc") {
    samples <- .ellipsoid_sampler(nsamples, boss$essential_support)
    log_unnormalized <- apply(samples, 1, log_density)
    unnormalized <- exp(log_unnormalized)
    avg <- mean(unnormalized)

    vol <- .ellipsoid_volume(boss$essential_support$H, boss$essential_support$chi2_radius)
    Z_hat <- avg * vol

    boss$normalized_posterior <- function(x) {
      exp(log_density(x)) / Z_hat
    }

  } else if (int_method == "aghq") {
    if (D > 3) {
      warning("AGHQ integration is not recommended for D > 3 due to exponential quadrature cost.")
    }

    if (is.null(boss$lower) || is.null(boss$upper)) {
      stop("For AGHQ, boss$lower and boss$upper must be defined.")
    }

    if (D == 1){
      bounds <- .extract_numeric_bounds_from_ellipsoid(boss$essential_support)
      lower <- bounds[1]
      upper <- bounds[2]
    }else{
      lower <- boss$lower
      upper <- boss$upper
    }

    lognormconst <- aghq_bounded(
      log_f = log_density,
      lower = lower,
      upper = upper,
      k = aghq_K,
      startingvalue = aghq_startingvalue,
      ...
    )

    boss$normalized_posterior <- function(x) {
      exp(log_density(x) - lognormconst)
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


#' @title Approximate log normalizing constant over bounded domain using AGHQ
#' @description
#' Approximates the log normalizing constant of a function defined on a bounded domain
#' [lower, upper]^D by applying a change-of-variable to the real line and using AGHQ.
#'
#' @param log_f Function taking a vector x ∈ [lower, upper]^D and returning log-density.
#' @param lower, upper Vectors of lower and upper bounds (length D).
#' @param k Number of quadrature points per dimension (e.g., 3–15).
#' @param startingvalue Optional starting value in transformed y-space (default is 0).
#' @return log normalizing constant (approximate)
#' @importFrom aghq aghq
#' @importFrom numDeriv grad hessian
#' @keywords internal
aghq_bounded <- function(log_f, lower, upper, k = 5, startingvalue = NULL) {
  stopifnot(length(lower) == length(upper))
  D <- length(lower)

  if (D > 3) {
    warning("aghq_bounded is only recommended for D ≤ 3 due to exponential growth of grid size.")
  }

  if (is.null(startingvalue)) {
    startingvalue <- rep(0, D)  # corresponds to midpoint in x-space
  }

  # transformed log-density in y-space
  transformed_log_f <- function(y) {
    y <- as.numeric(y)
    x <- pnorm(y) * (upper - lower) + lower
    log_jacobian <- sum(log(upper - lower)) + sum(dnorm(y, log = TRUE))
    log_f(x) + log_jacobian
  }

  ff <- list(
    fn = transformed_log_f,
    gr = function(y) numDeriv::grad(transformed_log_f, y),
    he = function(y) numDeriv::hessian(transformed_log_f, y)
  )

  aghq::aghq(ff = ff, k = k, startingvalue = startingvalue)$normalized_posterior$lognormconst
}





