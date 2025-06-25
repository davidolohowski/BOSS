#' Print Method for BOSS Objects
#'
#' Displays key information about a fitted BOSS object in a readable format.
#'
#' @param x A \code{boss} object.
#' @param ... Additional arguments (ignored).
#' @export
print.boss <- function(x, ...) {
  cat("Fitted BOSS Object\n")
  cat(strrep("-", 20), "\n")
  cat("Dimension of the conditioning parameter: ", x$D, "\n")
  cat("Total number of design points placed in the first stage: ",
      nrow(x$design_points$x_original), "\n")
  if (!is.null(x$essential_design_points$x_original)) {
    cat("Final number of design points after the second stage: ",
        nrow(x$essential_design_points$x_original), "\n")
  }
  fi <- if (!is.null(x$fill_in)) signif(x$fill_in, 6) else "NA"
  cat("Final fill-in distance:                 ", fi, "\n")
  mode_pt <- if (!is.null(x$mode)) {
    paste(signif(x$mode, 6), collapse = ", ")
  } else {
    "NA"
  }
  cat("Current mode of the conditioning param:  ", mode_pt, "\n")
  invisible(x)
}



#' Summary Method for BOSS Objects
#'
#' Provides a detailed summary of a fitted BOSS object, including:
#'   1) Initial design size,
#'   2) Fill-in stage details,
#'   3) Essential support properties,
#'   4) Inferred mode location.
#'
#' @param object A \code{boss} object.
#' @param ... Additional arguments (ignored).
#' @export
summary.boss <- function(object, ...) {
  cat("Detailed Summary of Fitted BOSS Object\n")
  cat(strrep("=", 36), "\n\n")

  ## 1) Initial design
  n_initial <- nrow(object$design_points$x_original)
  cat("1) Initial design:\n")
  cat("   - Number of design points placed in stage 1: ", n_initial, "\n\n")

  ## 2) Fill-in stage
  n_ess   <- if (!is.null(object$essential_design_points$x_original)) {
    nrow(object$essential_design_points$x_original)
  } else {
    0
  }
  n_added <- n_ess - n_initial
  fi      <- if (!is.null(object$fill_in)) signif(object$fill_in, 6) else NA
  cat("2) Fill-in stage:\n")
  cat("   - Points added in stage 2:           ", n_added, "\n")
  cat("   - Final fill-in distance achieved:   ", fi, "\n\n")

  ## 3) Essential support
  cat("3) Essential support:\n")
  cat("   - Total points in support:           ", n_ess, "\n")
  if (!is.null(object$essential_support)) {
    es    <- object$essential_support
    eigs  <- es$eig_values
    radii <- sqrt(es$chi2_radius / abs(eigs))
    minor <- signif(min(radii), 6)
    major <- signif(max(radii), 6)
    cat("   - Minor axis length:                ", minor, "\n")
    cat("   - Major axis length:                ", major, "\n\n")
  }

  ## 4) Inferred mode
  cat("4) Inferred mode:\n")
  if (!is.null(object$mode)) {
    mode_pt <- paste(signif(object$mode, 6), collapse = ", ")
    cat("   - Mode location in input space:      ", mode_pt, "\n")
  } else {
    cat("   - Mode location:                     NA\n")
  }

  invisible(object)
}





#' Plot Method for BOSS Objects with Optional Essential-Support Plot
#'
#' For D = 1 produces:
#' 1) Surrogate mean + design points + essential‐point ticks.
#' Then (without error) always draws:
#' 2) Essential-support schematic in 1D or 2D.
#'
#' @param x A \code{boss} object.
#' @param n_grid Integer; number of grid points for the surrogate (default 200).
#' @param tick_frac Numeric; fraction of the y-range for tick height (default 0.03).
#' @param ask Logical; if TRUE, waits for <Enter> before the second plot (default FALSE).
#' @param ... Additional arguments for the first \code{plot()} call.
#' @export
plot.boss <- function(x,
                      n_grid    = 200,
                      tick_frac = 0.03,
                      ask        = FALSE,
                      ...) {

  ## 1) If D==1, draw surrogate + design + ticks
  if (x$D == 1) {
    ys      <- x$essential_design_points$y
    ed_x    <- x$essential_design_points$x_original

    x_min <- min(ed_x); x_max <- max(ed_x)
    grid_x <- seq(x_min, x_max, length.out = n_grid)
    mu_hat <- sapply(grid_x, x$surrogate)

    plot(grid_x, mu_hat, type = "l", lwd = 2,
         xlab = "x", ylab = "Surrogate mean",
         main = "BOSS Surrogate and Essential Design Ticks", ...)
    points(ed_x, ys, pch = 16)
    usr   <- par("usr")
    y_bot <- usr[3]
    y_top <- usr[3] + tick_frac * (usr[4] - usr[3])
    segments(x0 = ed_x, y0 = y_bot, x1 = ed_x, y1 = y_top,
             col = "red", lwd = 2)
  }

  ## 2) Essential-support schematic
  if (ask) readline("Press <Enter> to see essential-support plot...")

  es     <- x$essential_support
  center <- es$center
  cval   <- es$chi2_radius

  if (length(center) == 1) {
    # 1D support interval
    H11    <- as.numeric(-es$H)
    radius <- if (H11 > 0) sqrt(cval / H11) else NA_real_
    left   <- center - radius
    right  <- center + radius

    plot(c(left, right), c(0, 0), type = "n", yaxt = "n",
         xlab = "x", ylab = "", main = "Essential Support & Points")
    segments(x0 = left, y0 = 0, x1 = right, y1 = 0, lwd = 3)
    points(center, 0, pch = 19)
    text(center, 0, "center", pos = 3)

    ed_x <- x$essential_design_points$x_original
    tick_h <- 0.1
    for (xi in ed_x) {
      segments(x0 = xi, y0 = -tick_h, x1 = xi, y1 = tick_h,
               col = "blue", lwd = 2)
    }

  } else if (length(center) == 2) {
    # 2D ellipse
    pts <- x$essential_design_points$x_original
    xlim <- range(pts[,1]); ylim <- range(pts[,2])
    plot(NA, xlim = xlim, ylim = ylim,
         xlab = "x1", ylab = "x2", main = "Essential Support & Points")
    points(pts, pch = 16)
    eig   <- eigen(es$H)
    vecs  <- eig$vectors; vals <- abs(eig$values)
    radii <- sqrt(cval / vals)
    angles <- seq(0, 2*pi, len = 200)
    unit   <- rbind(cos(angles), sin(angles))
    boundary <- t(vecs %*% diag(radii) %*% unit + center)
    lines(boundary, col = "red", lwd = 2)

  } else {
    warning("Essential-support plot only implemented for D=1 or D=2")
  }

  invisible(NULL)
}



