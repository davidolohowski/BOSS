#' Extract Interior and Boundary Design Points for a \code{boss} Object
#'
#' Given a \code{boss} object with \code{essential_support} computed (in original input space),
#' builds \code{essential_design_points}, including:
#'   1. all design points inside the ellipsoid,
#'   2. sampled points on the ellipsoid boundary,
#'   3. truncating any boundary samples to lie within the original \code{lower}–\code{upper} box.
#'
#' @param boss_result A \code{boss} object containing non-NULL
#'   \code{essential_support}, and with \code{lower} and \code{upper} fields.
#' @param boundary_initial Integer; number of boundary samples when \eqn{D >= 2} (default 100).
#'
#' @return The same \code{boss_result}, with
#'   \code{essential_design_points} list:
#'   \describe{
#'     \item{\code{x}}{Scaled points inside support + boundary samples (matrix).}
#'     \item{\code{x_original}}{Original-space points inside support + boundary samples (matrix).}
#'     \item{\code{y}}{Responses: real values for interior points, \code{NA} for boundary samples.}
#'   }
#'
#' @export
construct_essential_designs <- function(boss_result, boundary_initial = 100) {
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

  # 1. interior points
  diffs <- t(t(Xo) - mu)
  quad  <- rowSums((diffs %*% H) * diffs)
  idx   <- which(quad <= cval)
  Xs_int  <- Xs[idx,, drop=FALSE]
  Xo_int  <- Xo[idx,, drop=FALSE]
  Y_int   <- Y[idx]

  # 2. boundary samples
  if (D == 1) {
    H11 <- as.numeric(H)
    rad <- if (H11 > 0) sqrt(cval / H11) else NA_real_
    Xo_bdy <- matrix(mu + c(-rad, +rad), ncol = 1)
    Xs_bdy <- matrix(NA_real_, 2, 1)
    Y_bdy  <- rep(NA_real_, 2)
  } else {
    Vs <- matrix(rnorm(D * boundary_initial), ncol = D)
    Vs <- Vs / sqrt(rowSums(Vs^2))
    qs <- rowSums((Vs %*% H) * Vs)
    valid <- qs > 0
    Vs    <- Vs[valid,, drop=FALSE]
    ts    <- sqrt(cval / qs[valid])
    Xo_bdy <- sweep(Vs, 1, ts, "*") + matrix(mu, nrow = nrow(Vs), ncol = D, byrow = TRUE)
    Xs_bdy <- matrix(NA_real_, nrow(Xo_bdy), D)
    Y_bdy  <- rep(NA_real_, nrow(Xo_bdy))
  }

  # 3. combine
  Xo_all <- rbind(Xo_int, Xo_bdy)
  Xs_all <- rbind(Xs_int, Xs_bdy)
  Y_all  <- c(Y_int, Y_bdy)

  # 4. truncate to [lower, upper]
  lower <- boss_result$lower
  upper <- boss_result$upper
  # replicate bounds
  lb_mat <- matrix(lower, nrow = nrow(Xo_all), ncol = D, byrow = TRUE)
  ub_mat <- matrix(upper, nrow = nrow(Xo_all), ncol = D, byrow = TRUE)
  Xo_all <- pmin(pmax(Xo_all, lb_mat), ub_mat)

  # 5. recompute scaled coords
  Xs_all <- sweep(Xo_all, 2, lower, "-") / (upper - lower)

  boss_result$essential_design_points <- list(
    x          = Xs_all,
    x_original = Xo_all,
    y          = Y_all
  )
  return(boss_result)
}


#' Compute Fill-In Distance from Essential Design Points
#'
#' Given a \code{boss_result} that already has
#' \code{essential_design_points} (from \code{construct_essential_designs}),
#' compute the maximal minimal Euclidean distance among those points.
#'
#' @param boss_result A \code{boss} object containing a non-NULL
#'   \code{essential_design_points} list with elements \code{x}, \code{x_original}, and \code{y}.
#'
#' @return The input \code{boss_result}, augmented with a numeric
#'   component \code{fill_in}, the maximum of each point’s nearest-neighbor distance.
#'
#' @importFrom stats dist
#' @export
compute_fill_in <- function(boss_result) {
  ed <- boss_result$essential_design_points
  if (is.null(ed)) {
    stop("essential_design_points must be constructed first")
  }

  X_pool <- ed$x_original
  n_pool <- nrow(X_pool)
  if (n_pool < 2) {
    boss_result$fill_in <- 0
    return(boss_result)
  }

  # all pairwise distances
  dmat <- as.matrix(stats::dist(X_pool))
  diag(dmat) <- Inf
  min_d <- apply(dmat, 1, min, na.rm = TRUE)

  boss_result$fill_in <- max(min_d, na.rm = TRUE)
  return(boss_result)
}


#' Greedy Fill-In of Essential Design Points to Reach Target Spacing (with Verbosity)
#'
#' Starting from a \code{boss_result} that already has
#' \code{essential_design_points$x_original} and a computed \code{fill_in},
#' iteratively add new points at midpoints of the largest gaps until
#' \code{fill_in <= h} or until \code{max_add} points have been added.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_design_points$x_original} and numeric \code{fill_in}.
#' @param h Numeric; target fill-in distance (must be \(\ge0\)).
#' @param max_add Integer; maximum number of new points to add (default 20).
#' @param verbose Integer (0, 1, or 2); level of printed feedback:
#'   \describe{
#'     \item{0}{Silent.}
#'     \item{1}{Print total iterations and final fill_in.}
#'     \item{2}{Also print fill_in after each addition.}
#'   }
#'
#' @return The same \code{boss_result}, with its
#'   \code{essential_design_points$x_original} augmented by new points,
#'   and its \code{fill_in} updated to \(\le h\) if possible.
#'
#' @importFrom stats dist
#' @export
fill_in <- function(boss_result, h, max_add = 20, verbose = 0) {
  ed <- boss_result$essential_design_points
  if (is.null(ed) || is.null(ed$x_original)) {
    stop("essential_design_points$x_original must exist")
  }
  if (!is.numeric(boss_result$fill_in)) {
    stop("Please compute fill_in (via compute_fill_in) before calling fill_in()")
  }
  if (!is.numeric(h) || length(h) != 1 || h < 0) {
    stop("`h` must be a non-negative numeric scalar")
  }
  if (!verbose %in% 0:2) {
    stop("`verbose` must be 0, 1, or 2")
  }

  Xo    <- ed$x_original
  Xs    <- ed$x
  y     <- ed$y
  added <- 0

  while (boss_result$fill_in > h && added < max_add) {
    # compute distances
    dmat    <- as.matrix(stats::dist(Xo))
    diag(dmat) <- Inf
    nn_dist <- apply(dmat, 1, min, na.rm = TRUE)

    # find largest gap
    i_max <- which.max(nn_dist)
    j_max <- which.min(dmat[i_max, ])

    # midpoint
    new_pt <- as.numeric((Xo[i_max, ] + Xo[j_max, ]) / 2)

    # append
    Xo <- rbind(Xo, matrix(new_pt, nrow = 1))
    rownames(Xo) <- NULL
    Xs <- rbind(Xs, (new_pt - boss_result$lower) / (boss_result$upper - boss_result$lower))
    y <- c(y, NA_real_)  # new point has no response yet
    added <- added + 1

    # update and recompute fill_in
    boss_result$essential_design_points$x_original <- Xo
    boss_result$essential_design_points$x <- Xs
    boss_result$essential_design_points$y <- y
    boss_result <- compute_fill_in(boss_result)

    if (verbose == 2) {
      cat("Iteration", added, ": new fill_in =", signif(boss_result$fill_in, 6), "\n")
    }
  }

  if (verbose >= 1) {
    cat("fill_in completed:", added, "point(s) added; final fill_in =",
        signif(boss_result$fill_in, 6), "\n")
  }
  if (boss_result$fill_in > h) {
    warning("Reached max_add without achieving target fill-in")
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









