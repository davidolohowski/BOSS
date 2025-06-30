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


#' Approximate Fill-In Distance via Lattice Interior + Uniform Boundary Sampling (with Truncation)
#'
#' Estimate the fill-in distance
#' \(\sup_{x\in\Omega}\min_j \|x - x_j\|\)
#' by:
#' 1. Constructing a regular lattice inside the ellipsoid Ω,
#' 2. Sampling points uniformly on the ellipsoid boundary,
#' 3. Truncating all candidates to lie within the original \code{lower}/\code{upper} box,
#' 4. Computing for each candidate its nearest-design-point distance,
#' 5. Taking the maximum of these minima.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support}, \code{essential_design_points\$x_original},
#'   and containing \code{lower} and \code{upper}.
#' @param grid_n Integer; number of equally spaced points per axis for the interior lattice (default 200).
#' @param boundary_n Integer; number of random samples on the boundary (default 200).
#'
#' @return The input \code{boss_result}, with \code{fill_in} updated.
#'
#' @importFrom stats dist
#' @export
compute_fill_in <- function(boss_result, grid_n = 200, boundary_n = 200) {
  ess <- boss_result$essential_support
  ed  <- boss_result$essential_design_points
  if (is.null(ess) || is.null(ed$x_original)) {
    stop("Please compute essential_support and construct_essential_designs first.")
  }
  # unpack
  center <- ess$center
  Hinv   <- -ess$H            # flipped to PSD
  cval   <- ess$chi2_radius
  Xd     <- ed$x_original
  lower  <- boss_result$lower
  upper  <- boss_result$upper
  D      <- length(center)

  # 1) axis-aligned bounding box of Ω
  eig      <- eigen(Hinv, symmetric = TRUE)
  vecs     <- eig$vectors
  vals     <- eig$values
  radii    <- sqrt(cval / vals)
  abs_vecs <- abs(vecs %*% diag(radii))
  half_box <- rowSums(abs_vecs)
  mins     <- center - half_box
  maxs     <- center + half_box

  # 2) build interior lattice and keep only points in Ω
  seqs   <- lapply(seq_len(D), function(i)
    seq(mins[i], maxs[i], length.out = grid_n))
  grid   <- as.matrix(do.call(expand.grid, seqs))
  diffs  <- sweep(grid, 2, center, "-")
  in_idx <- rowSums((diffs %*% Hinv) * diffs) <= cval
  Xi     <- grid[in_idx, , drop = FALSE]

  # 3) uniform boundary sampling
  if (D == 1) {
    H11 <- as.numeric(Hinv)
    r   <- if (H11 > 0) sqrt(cval / H11) else NA_real_
    Xb  <- matrix(center + c(-r, +r), ncol = 1)
  } else {
    Ainv   <- solve(Hinv)
    W      <- MASS::mvrnorm(boundary_n, mu = rep(0, D), Sigma = Ainv)
    norms2 <- rowSums((W %*% Hinv) * W)
    Xb     <- sweep(W, 1, sqrt(cval / norms2), "*")
    Xb     <- sweep(Xb, 2, center, "+")
  }

  # 4) combine and truncate to [lower, upper]
  Xc <- rbind(Xi, Xb)
  lb_mat <- matrix(lower, nrow = nrow(Xc), ncol = D, byrow = TRUE)
  ub_mat <- matrix(upper, nrow = nrow(Xc), ncol = D, byrow = TRUE)
  Xc <- pmin(pmax(Xc, lb_mat), ub_mat)

  # 5) compute minimal distance of each candidate to existing design
  dmat <- as.matrix(stats::dist(rbind(Xd, Xc)))
  K    <- nrow(Xd)
  dsub <- dmat[(K+1):nrow(dmat), 1:K, drop = FALSE]
  min_d <- apply(dsub, 1, min)

  # 6) fill_in = maximum of these minimal distances
  boss_result$fill_in <- max(min_d, na.rm = TRUE)
  return(boss_result)
}



#' Fill-In by Greedy Selection with Candidate Sufficiency Checks
#'
#' At each iteration:
#' 1. Builds a fixed pool of candidate points (interior lattice + uniform boundary),
#'    truncated to the original box.
#' 2. Checks that the candidate pool size is at least \code{max_add}.
#' 3. Greedily selects the candidate farthest from the current design,
#'    appends it, and recomputes \code{fill_in}.
#' 4. Stops when \code{fill_in <= h} or \code{max_add} points have been added.
#'
#' If the candidate pool is too small to potentially reach \code{max_add},
#' or if even the best candidate cannot reduce \code{fill_in} below \code{h},
#' a warning is issued.
#'
#' @param boss_result A \code{boss} object with non-NULL
#'   \code{essential_support} and \code{essential_design_points\$x_original},
#'   and containing \code{lower} and \code{upper}.
#' @param h Numeric; target fill-in distance (>= 0).
#' @param max_add Integer; maximum number of points to add (default 100).
#' @param grid_n Integer; resolution per axis for interior lattice (default 200).
#' @param boundary_n Integer; number of boundary samples (default 200).
#' @param verbose Integer in \{0,1,2\}; 0=silent, 1=summary, 2=per-iteration.
#'
#' @return The updated \code{boss_result}, with its
#'   \code{essential_design_points} and \code{fill_in} refreshed.
#'
#' @importFrom stats dist
#' @export
fill_in <- function(boss_result,
                    h          ,
                    max_add    = 100,
                    grid_n     = 200,
                    boundary_n = 200,
                    verbose    = 0) {
  # prerequisites
  ess <- boss_result$essential_support
  ed  <- boss_result$essential_design_points
  if (is.null(ess) || is.null(ed$x_original)) {
    stop("Please run compute_essential_support() and construct_essential_designs() first")
  }
  if (!is.numeric(boss_result$fill_in)) {
    stop("Please compute fill_in() before calling fill_in()")
  }
  if (!is.numeric(h) || length(h)!=1 || h<0) {
    stop("`h` must be a non-negative numeric scalar")
  }
  if (!verbose %in% 0:2) {
    stop("`verbose` must be 0, 1, or 2")
  }

  # unpack
  center <- ess$center
  Hinv   <- -ess$H
  cval   <- ess$chi2_radius
  Xo     <- ed$x_original
  Xs     <- ed$x
  Y      <- ed$y
  lower  <- boss_result$lower
  upper  <- boss_result$upper
  D      <- length(center)
  added  <- 0

  # 1) Build interior lattice
  eig   <- eigen(Hinv, symmetric=TRUE)
  vecs  <- eig$vectors; vals <- eig$values
  radii <- sqrt(cval / vals)
  halfB <- rowSums(abs(vecs %*% diag(radii)))
  mins  <- center - halfB
  maxs  <- center + halfB
  seqs  <- lapply(seq_len(D),
                  function(i) seq(mins[i], maxs[i], length.out = grid_n))
  grid  <- as.matrix(do.call(expand.grid, seqs))
  diffs <- sweep(grid, 2, center, "-")
  Xi    <- grid[rowSums((diffs %*% Hinv) * diffs) <= cval, , drop=FALSE]

  # 2) Uniform boundary sampling (area-uniform)
  if (D == 1) {
    H11 <- as.numeric(Hinv)
    r   <- if (H11 > 0) sqrt(cval / H11) else NA_real_
    Xb  <- matrix(center + c(-r, r), ncol=1)
  } else {
    Ainv   <- solve(Hinv)
    W      <- MASS::mvrnorm(boundary_n, mu=rep(0,D), Sigma=Ainv)
    norms2 <- rowSums((W %*% Hinv)*W)
    Xb     <- sweep(W, 1, sqrt(cval / norms2), "*")
    Xb     <- sweep(Xb, 2, center, "+")
  }

  # 3) Pool and truncate to [lower, upper]
  candidates <- unique(
    pmin(pmax(rbind(Xi, Xb),
              matrix(lower, nrow=nrow(rbind(Xi,Xb)), ncol=D, byrow=TRUE)),
         matrix(upper, nrow=nrow(rbind(Xi,Xb)), ncol=D, byrow=TRUE))
  )

  # 4) Check candidate pool size
  n_cand <- nrow(candidates)
  if (n_cand < max_add) {
    warning(sprintf("Candidate pool size (%d) < max_add (%d); cannot add that many points, consider increasing grid_n and boundary_n.",
                    n_cand, max_add))
  }


  dc_full <- as.matrix(stats::dist(unique(rbind(Xo, candidates))))
  dc_full[dc_full == 0] <- Inf  # avoid zero distances
  min_d_achievable <- apply(dc_full, 1, min)
  best_achievable_dist <- max(min_d_achievable)
  if (best_achievable_dist > h) {
    warning(
      sprintf(
        "Fill-in distance between candidate points is larger than the target h! \n  (%.6f) > target h (%.6f); cannot reach h, consider increasing grid_n and boundary_n.",
        best_achievable_dist,
        h
      )
    )
  }

  # 5) Greedy addition
  while (boss_result$fill_in > h && added < max_add && nrow(candidates) > 0) {
    dmat_cand <- as.matrix(stats::dist(rbind(Xo, candidates)))
    K <- nrow(Xo)
    dc <- dmat_cand[(K+1):nrow(dmat_cand), 1:K, drop=FALSE]
    min_d <- apply(dc, 1, min)
    best_idx <- which.max(min_d)
    best_dist <- min_d[best_idx]

    # append best
    new_o <- candidates[best_idx, ]
    new_s <- (new_o - lower) / (upper - lower)
    Xo   <- rbind(Xo, new_o)
    Xs   <- rbind(Xs, new_s)
    Y    <- c(Y, NA_real_)
    added <- added + 1

    # remove used candidate
    candidates <- candidates[-best_idx, , drop=FALSE]

    # update
    boss_result$essential_design_points$x_original <- Xo
    boss_result$essential_design_points$x          <- Xs
    boss_result$essential_design_points$y          <- Y

    # multiply by 10 to increase resolution
    boss_result <- compute_fill_in(boss_result,
                                   grid_n = (grid_n*10), boundary_n = (boundary_n*10))

    if (verbose == 2) {
      cat("Iteration", added,
          ": added point =", signif(new_o, 6),
          "-> fill_in =", signif(boss_result$fill_in, 6), "\n")
    }
  }

  # 6) Final message
  if (verbose >= 1) {
    cat("fill_in completed: added", added,
        "point(s); final fill_in =", signif(boss_result$fill_in, 6), "\n")
  }
  if (boss_result$fill_in > h) {
    warning("Reached termination without achieving target fill-in")
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









