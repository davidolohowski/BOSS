#' Compute Pairwise Squared Distances
#'
#' This function computes the matrix of squared Euclidean distances
#' between the rows of two numeric inputs.  If either \code{X} or \code{Y}
#' is a vector, it will be treated as a single-column matrix.
#'
#' @param X A numeric matrix or vector.  Rows represent observations.
#' @param Y A numeric matrix or vector.  Rows represent observations.
#' @return A numeric matrix of size \code{nrow(X) x nrow(Y)}, where entry
#'   \code{[i, j]} is the squared distance between row \code{i} of \code{X}
#'   and row \code{j} of \code{Y}.
#' @details Internally, uses the identity
#'   \\[||x - y||^2 = ||x||^2 + ||y||^2 - 2 x^\top y\\]
#'   to compute all pairwise squared distances efficiently.
#' @examples
#' # Two 2D points (Not run)
#' #X <- matrix(c(0,0, 1,1), nrow = 2, byrow = TRUE)
#' # Y <- matrix(c(1,0, 0,1), nrow = 2, byrow = TRUE)
#' # compute_sq_dist(X, Y)
#'
#' # Compare to manual calculation
#' # (0,0) vs (1,0): 1^2 + 0^2 = 1
#' # (0,0) vs (0,1): 0^2 + 1^2 = 1
compute_sq_dist <- function(X, Y) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  XY <- tcrossprod(X, Y)
  sq_dist <- matrix(XX, ncol = nrow(Y), nrow = nrow(X)) +
    t(matrix(YY, ncol = nrow(X), nrow = nrow(Y))) -
    2 * XY
  return(sq_dist)
}




#' Generate a General SE/Matérn Covariance Function
#'
#' Creates a covariance function that computes either the Squared Exponential (SE)
#' kernel (when \code{nu = Inf}) or the Matérn kernel (for finite \code{nu})
#' between two sets of input points in \eqn{d}-dimensional space.
#'
#' @param length_scale A positive numeric scalar or vector.
#'   If scalar, the same length-scale is used for all dimensions.
#'   If a vector, its length must equal the number of columns of the input data.
#' @param signal_var Positive numeric scalar specifying the signal variance (\eqn{\sigma^2}).
#' @param nu Positive numeric specifying the smoothness parameter for the Matérn kernel.
#'   Use \code{Inf} to select the Squared Exponential kernel.
#'
#' @return A function \code{cov_fun(x, x_prime)} which accepts:
#'   \itemize{
#'     \item \code{x}: numeric matrix or vector (observations in rows).
#'     \item \code{x_prime}: numeric matrix or vector (observations in rows).
#'   }
#'   and returns the covariance matrix between \code{x} and \code{x_prime}.
#'
#' @details
#' The returned function first rescales each dimension by \code{length_scale},
#' then computes the pairwise squared Euclidean distances.
#'
#' - If \code{nu = Inf}, it returns the SE kernel:
#'   \deqn{K_{ij} = \sigma^2 \exp\bigl(-\tfrac{1}{2} \|x_i - x'_j\|^2\bigr).}
#' - Otherwise, it constructs the Matérn kernel with smoothness \code{nu}:
#'   \deqn{K_{ij} = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)} \bigl(\sqrt{2\nu}\,r_{ij}\bigr)^\nu
#'     K_\nu\bigl(\sqrt{2\nu}\,r_{ij}\bigr),}
#'   where \eqn{r_{ij} = \|x_i - x'_j\|}.
#' Internally, for half-integer \code{nu} this is computed via an exact polynomial–exponential form.
#'
#' @examples
#' # 1D SE kernel (Not run)
#' # se_cov <- cov_generator(length_scale = 0.5, signal_var = 2, nu = Inf)
#' # x <- seq(0, 1, length.out = 5)
#' # se_cov(x, x)
#'
#' # 2D Matérn kernel with nu = 1.5
#' # matern_cov <- cov_generator(length_scale = c(1, 2), signal_var = 1, nu = 1.5)
#' # X <- matrix(runif(10), ncol = 2)
#' # matern_cov(X, X)
#'
#' @export
cov_generator <- function(length_scale = 1, signal_var = 1, nu = Inf) {

  cov_fun <- function(x, x_prime) {
    # ensure matrix form
    if (is.vector(x))       x       <- matrix(x,       ncol = 1)
    if (is.vector(x_prime)) x_prime <- matrix(x_prime, ncol = 1)
    if (ncol(x) != ncol(x_prime))
      stop("x and x_prime must have same number of columns")

    d <- ncol(x)
    # build a length-scale vector
    if (length(length_scale) == 1) {
      ls <- rep(length_scale, d)
    } else {
      ls <- length_scale
    }
    if (length(ls) != d)
      stop("length_scale must be length 1 or match ncol(x)")

    # pre-scale each column by its length-scale
    x_s       <- sweep(x,       2, ls, "/")
    x_prime_s <- sweep(x_prime, 2, ls, "/")

    # compute scaled Euclidean distances
    sq_dist_s <- compute_sq_dist(x_s, x_prime_s)

    # if same input, mirror to keep it symmetric
    if (identical(x, x_prime)) {
      ut <- upper.tri(sq_dist_s)
      sq_dist_s[ut] <- t(sq_dist_s)[ut]
    }

    # Safe correction for numerical errors
    sq_dist_s[sq_dist_s < 0] <- 0

    # Decide based on nu
    if (is.infinite(nu)) {
      # Square Exponential
      cov_matrix <- signal_var * exp(-sq_dist_s / 2)

    } else {
      # Matérn
      m <- ceiling(nu) - 1
      nu_adj <- m + 0.5

      dist_s <- sqrt(sq_dist_s)
      sqrt_2nu <- sqrt(2 * nu_adj)
      z <- sqrt_2nu * dist_s
      exp_part <- exp(-z)

      k_vec <- 0:m
      coeff_vec <- exp(lgamma(m + 1) - lgamma(2*m + 1) +
                         lgamma(m + k_vec + 1) -
                         (lgamma(k_vec + 1) + lgamma(m - k_vec + 1)))

      # Build polynomial part
      z_vec <- as.vector(2*z)
      z_powers <- outer(z_vec, m - k_vec, "^")
      poly_part_vec <- z_powers %*% coeff_vec
      poly_part <- matrix(poly_part_vec, nrow = nrow(z), ncol = ncol(z))

      cov_matrix <- signal_var * poly_part * exp_part
    }

    return(t(cov_matrix))
  }

  return(cov_fun)
}



#' Gaussian Process Prediction (with Optional Quadratic Mean)
#'
#' Compute the posterior mean, covariance, and a single simulation draw
#' from a Gaussian Process (GP) posterior at new input locations.
#' Supports both zero-mean GP and a linear-plus-quadratic mean model.
#'
#' @param data A list with elements:
#'   \itemize{
#'     \item \code{x}: numeric vector or matrix of observed inputs (rows = observations).
#'     \item \code{y}: numeric vector of observed responses.
#'   }
#' @param x_pred Numeric vector or matrix of input locations at which to predict.
#'   Must have the same number of columns as \code{data$x} if \code{data$x} is a matrix.
#' @param noise_var Positive numeric scalar.  Observation noise variance to add
#'   on the diagonal of the training covariance matrix.
#' @param choice_cov A covariance-generating function, such as one returned by
#'   \code{\link[BOSS]{cov_generator}}, that accepts two inputs
#'   \code{(x1, x2)} and returns their covariance matrix.
#' @param quad Logical; if \code{FALSE} (default), use a zero-mean GP model.
#'   If \code{TRUE}, fit a linear-plus-quadratic mean model
#'   \eqn{\beta_0 + x^\top \beta_1 + (x \otimes x)^\top \beta_2} alongside the GP.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{x}}{The prediction inputs \code{x_pred}.}
#'     \item{\code{mean}}{Numeric vector of posterior mean values at \code{x_pred}.}
#'     \item{\code{var}}{Posterior covariance matrix among \code{x_pred}.}
#'     \item{\code{sim}}{A single draw from the multivariate normal with
#'       mean \code{mean} and covariance \code{var}.}
#'   }
#'
#' @details
#' For the zero-mean case (\code{quad = FALSE}), the function:
#' \itemize{
#'   \item Constructs \eqn{K_{oo}, K_{op}, K_{pp}} via \code{choice_cov}.
#'   \item Adds \code{noise_var} to the diagonal of \eqn{K_{oo}} and computes
#'     its Cholesky factor \code{L}.
#'   \item Computes the conditional mean \eqn{K_{op} L^{-T} L^{-1} y} and
#'     covariance \eqn{K_{pp} - K_{op} L^{-T} L^{-1} K_{op}^T}.
#' }
#' If \code{quad = TRUE}, a design matrix with intercept, linear, and
#' unique quadratic terms is built for both observed and prediction inputs,
#' and coefficients \eqn{\beta} are estimated via generalized least squares.
#' The posterior mean and covariance combine both the mean model and GP residuals.
#'
#' @examples
#' # Zero-mean GP with SE kernel (Not run)
#' # cov_fn <- cov_generator(length_scale = 1, signal_var = 1, nu = Inf)
#' # data <- list(x = seq(0, 1, length.out = 10),
#' #               y = sin(pi * seq(0, 1, length.out = 10)))
#' # preds <- predict_gp(data, x_pred = seq(0, 1, length.out = 50),
#' #                      noise_var = 1e-4, choice_cov = cov_fn)
#' # plot(preds$x, preds$mean, type = "l")
#' #
#' # GP with quadratic mean
#' # preds_quad <- predict_gp(data, x_pred = seq(0, 1, length.out = 50),
#' #                         noise_var = 1e-4, choice_cov = cov_fn,
#' #                       quad = TRUE)
#'
#' @importFrom MASS mvrnorm
#' @export
predict_gp <- function(data,
                       x_pred,
                       noise_var = 1e-6,
                       choice_cov,
                       quad = FALSE) {
  # Extract x and y
  x_obs <- data$x
  y_obs <- data$y
  if (is.vector(x_obs)) {
    N <- length(x_obs); D <- 1
  } else {
    N <- nrow(x_obs);   D <- ncol(x_obs)
  }

  if (!quad) {
    # Zero-mean GP prediction
    K_obs_obs   <- choice_cov(x_obs,  x_obs)
    K_obs_pred  <- choice_cov(x_obs,  x_pred)
    K_pred_pred <- choice_cov(x_pred, x_pred)

    K_obs_obs <- K_obs_obs + noise_var * diag(N)
    L         <- chol(K_obs_obs)

    # Posterior mean
    Ly        <- forwardsolve(t(L), y_obs)
    cond_mean <- K_obs_pred %*% backsolve(L, Ly)

    # Posterior covariance
    LK        <- forwardsolve(t(L), t(K_obs_pred))
    cond_var  <- K_pred_pred - crossprod(LK)

    # Single sample
    sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)
    return(list(
      x    = x_pred,
      mean = as.vector(cond_mean),
      var  = cond_var,
      sim  = sim
    ))
  } else {
    # GP with quadratic mean
    K_obs_obs   <- choice_cov(x_obs, x_obs)
    K_obs_pred  <- choice_cov(x_obs, x_pred)
    K_pred_pred <- choice_cov(x_pred, x_pred)

    K_obs_obs <- K_obs_obs + noise_var * diag(N)
    L <- chol(K_obs_obs)

    # Design matrices: intercept + linear + unique quadratic terms
    X_covariate <- cbind(1, x_obs,
                         t(apply(x_obs, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))
    X_star <- cbind(1, x_pred,
                    t(apply(x_pred, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))

    # Estimate beta via GLS
    LXs <- forwardsolve(t(L), X_covariate)
    A   <- crossprod(LXs)

    tmp  <- forwardsolve(t(L), y_obs)
    Ly   <- backsolve(L, tmp)

    B   <- crossprod(X_covariate, Ly)
    beta <- solve(A, B)

    # Residual contributions
    residual  <- y_obs - X_covariate %*% beta
    tmp <- forwardsolve(t(L), residual)
    Lres <- backsolve(L, tmp)

    cond_mean <- X_star %*% beta + K_obs_pred %*% Lres

    # Covariance decomposition
    LK_pred   <- forwardsolve(t(L), t(K_obs_pred))
    var_part1 <- K_pred_pred - crossprod(LK_pred)
    delta     <- X_star - t(crossprod(LXs, LK_pred))
    chol_A    <- chol(A)
    Ldelta    <- forwardsolve(t(chol_A), t(delta))
    var_part2 <- crossprod(Ldelta)
    cond_var  <- var_part1 + var_part2

    # Single sample
    sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)

    return(list(
      x    = x_pred,
      mean = as.vector(cond_mean),
      var  = cond_var,
      sim  = sim
    ))
  }
}



#' Compute Gaussian Process Log-Likelihood (with Optional Quadratic Mean)
#'
#' Evaluate the (negative) log-likelihood of observed data under a Gaussian Process
#' model, optionally including a linear-plus-quadratic mean function.
#'
#' @param length_scale Positive numeric scalar or vector of length-scales for each dimension.
#' @param x Numeric vector or matrix of training inputs (rows = observations).
#' @param y Numeric vector of training responses.
#' @param signal_var Positive numeric scalar specifying the GP signal variance (\eqn{\sigma^2}).
#' @param noise_var Positive numeric scalar specifying observation noise variance.
#' @param D Integer. Dimensionality of the input space (number of columns in \code{x}).
#' @param quad Logical; if \code{FALSE} (default), assume zero-mean GP. If \code{TRUE},
#'   include a mean function \eqn{\beta_0 + x^\top \beta_1 + (x \otimes x)^\top \beta_2}.
#' @param nu Positive numeric specifying the smoothness parameter for the Matérn kernel.
#'   Defaults to \code{Inf}, selecting the Squared Exponential kernel.
#'
#' @return A numeric scalar giving the negative log-likelihood (up to an additive constant).
#'   Extremely large values (\eqn{>1e^{20}}) are returned if the computation produces NaN, NA,
#'   or infinite results.
#'
#' @details
#' **Zero-mean GP** (\code{quad = FALSE}):
#'
#' * Build covariance \eqn{C = K_{xx} + \sigma^2_n I\} via
#'   \code{\link[BOSS]{cov_generator}}.
#'
#' * Compute Cholesky factor \eqn{L L^T = C}.
#'
#' * Log-determinant: \eqn{\log\det(C) = 2\sum\log(\mathrm{diag}(L))}.
#' * Precision \eqn{Q = C^{-1} = \texttt{chol2inv}(L)}.
#'
#' * Negative log-likelihood:
#'   \deqn{
#'     \frac{1}{2}\,y^T Q\,y \;-\; \log\det(Q)
#'     \;-\;\sum \log\bigl[\text{LogNormal}(length\_scale\mid prior\_l)\bigr].
#'   }
#'
#' **GP with quadratic mean** (\code{quad = TRUE}):
#'
#' * Construct design matrix \eqn{X = [1,\,x,\,\mathrm{unique}(x \otimes x)]}.
#'
#' * Build \eqn{C = K_{xx} + \sigma^2_n I} and its precision \eqn{Q}.
#'
#' * Form marginal precision \eqn{S = X^\top Q X} and invert via Cholesky.
#'
#' * Compute two quadratic terms:
#'   \eqn{\tfrac12\,y^\top Q\,y} and \eqn{\tfrac12\,y^\top Q X S^{-1} X^\top Q\,y}.
#'
#' * Combine with log-determinants of \eqn{Q} and \eqn{S^{-1}} and the log-prior
#'   to yield the final negative log-likelihood.
#'
#' @examples
#' # Zero-mean SE kernel likelihood (Not run)
#' # x <- matrix(seq(0,1,length.out=5), ncol=1)
#' # y <- sin(2*pi*x)
#' # like0 <- compute_like(length_scale = 0.2,
#' #                        x = x, y = y,
#' #                        signal_var = 1,
#' #                        noise_var  = 1e-3,
#' #                        D = 1,
#' #                        prior_l_mean = 0,
#' #                        prior_l_sd   = 1,
#' #                        quad = FALSE,
#' #                        nu   = Inf)
#'
#' # GP with quadratic mean
#' # like2 <- compute_like(length_scale = 0.3,
#' #                       x = x, y = y,
#' #                        signal_var = 1,
#' #                        noise_var  = 1e-3,
#' #                        D = 1,
#' #                        prior_l_mean = 0,
#' #                        prior_l_sd   = 1,
#' #                        quad = TRUE,
#' #                        nu   = Inf)
#'
#' @export
compute_like <- function(length_scale, x, y,
                         signal_var, noise_var,
                         D,
                         quad = FALSE, nu = Inf) {
  if (!quad) {
    choice_cov <- cov_generator(length_scale = length_scale,
                                signal_var   = signal_var,
                                nu           = nu)
    C <- choice_cov(x, x)
    A <- C + noise_var * diag(nrow(x))
    L <- chol(A)

    # log-det and precision
    log_det_A <- sum(log(diag(L)))
    log_det_Q <- -log_det_A
    Q         <- chol2inv(L)

    like <- as.numeric((t(y) %*% Q %*% y) / 2
                       - log_det_Q)
    if (is.nan(like) || is.na(like)) return(1e20)
    if (like == -Inf)               return(-1e20)
    if (like ==  Inf)               return(1e20)
    return(like)

  } else {
    # design matrix with intercept, linear, and unique quadratic terms
    X_covariate <- cbind(
      rep(1, nrow(x)),
      x,
      t(apply(x, 1, function(y) ((y %*% t(y))[upper.tri(diag(ncol(x)), diag = TRUE)])))
    )

    choice_cov <- cov_generator(length_scale = length_scale,
                                signal_var   = signal_var,
                                nu           = nu)
    C <- choice_cov(x, x) + noise_var * diag(nrow(x))

    # precision and marginal precision S = X^T Q X
    Lc  <- chol(C)
    Q   <- chol2inv(Lc)
    S_mat <- crossprod(X_covariate, Q %*% X_covariate)
    Ls    <- chol(S_mat)
    KK    <- chol2inv(Ls)

    # log-determinants
    log_det_Q  <- -sum(log(diag(Lc)))
    log_det_KK <- -sum(log(diag(Ls)))

    # quadratic forms
    q1 <- as.numeric((t(y) %*% Q %*% y) / 2)
    q2 <- as.numeric((t(y) %*% Q %*% X_covariate %*% KK %*%
                        t(X_covariate) %*% Q %*% y) / 2)

    like <- q1 - q2 - log_det_Q + log_det_KK

    if (is.nan(like) || is.na(like)) return(1e20)
    if (like == -Inf)               return(-1e20)
    if (like ==  Inf)               return(1e20)
    return(like)
  }
}




#' Upper Confidence Bound (UCB) Acquisition Function
#'
#' Compute the UCB acquisition value (to be maximized) at candidate inputs for Bayesian optimization.
#'
#' @param x Numeric vector or matrix of candidate input locations.
#' @param data A list with observed GP data, containing elements \code{x} (inputs) and \code{y} (responses).
#' @param cov A covariance-generating function (e.g., returned by \code{\link[BOSS]{cov_generator}}).
#' @param nv Numeric scalar. Observation noise variance.
#' @param D Integer. Total number of GP training points.
#' @param d Integer. Dimensionality of the input space.
#' @param quad Logical; passed to \code{\link[BOSS]{predict_gp}}. If \code{TRUE}, uses a quadratic mean model.
#'
#' @return Numeric vector of UCB values (higher is better). Internally returns \(-\mu(x) - \sqrt{\beta\,\Sigma(x)}\)
#'   so that minimizing this corresponds to maximizing the classic UCB acquisition \(\mu + \sqrt{\beta\,\sigma^2}\).
#'
#' @details
#' The exploration–exploitation tradeoff parameter \(\beta\) is set as
#' \deqn{\beta = 2\log\bigl(D^2\pi^2/(6\,d)\bigr).}
#' The function calls \code{\link[BOSS]{predict_gp}} to obtain the posterior mean \(\mu(x)\)
#' and variance \(\sigma^2(x)\), then computes
#' \[
#'   \mathrm{UCB}(x) \;=\; \mu(x) \;+\;\sqrt{\beta\,\sigma^2(x)}.
#' \]
#'
#' @examples
#' # Not run
#' # cov_fn <- cov_generator(length_scale = 1, signal_var = 1, nu = Inf)
#' # data <- list(x = matrix(runif(20), ncol = 2),
#' #               y = rnorm(10))
#' # UCB(seq(0,1,length.out=5), data, cov_fn, nv = 1e-3, D = 10, d = 2)
#'
#' @export
UCB <- function(x, data, cov, nv, D, d, quad = FALSE){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv, quad = quad)

  # Compute the UCB acquisition function
  beta <- 2*log((D^2)*(pi^2)/(6*d))
  return(as.numeric(-fnew$mean - sqrt(beta*fnew$var)))
}




#' Pure Exploration Acquisition Function
#'
#' Compute a variance-only acquisition value for Bayesian optimization,
#' favoring points with the highest posterior uncertainty.
#'
#' @param x Numeric vector or matrix of candidate input locations.
#' @param data A list with observed GP data, containing elements \code{x} (inputs) and \code{y} (responses).
#' @param cov A covariance-generating function (e.g., returned by \code{\link[BOSS]{cov_generator}}).
#' @param nv Numeric scalar. Observation noise variance.
#' @param quad Logical; passed to \code{\link[BOSS]{predict_gp}}. If \code{TRUE}, uses a quadratic mean model.
#'
#' @return Numeric vector of exploration scores \(-\sigma^2(x)\). Minimizing this corresponds to selecting
#'   inputs with the largest posterior variance.
#'
#' @examples
#' # Not run
#' # cov_fn <- cov_generator(length_scale = 1, signal_var = 1, nu = Inf)
#' # data <- list(x = matrix(runif(20), ncol = 2),
#' #             y = rnorm(10))
#' # EXPLORE(seq(0,1,length.out=5), data, cov_fn, nv = 1e-3)
#'
#' @export
EXPLORE <- function(x, data, cov, nv, quad = FALSE){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv, quad = quad)

  # Compute the pure exploration acquisition function
  return(as.numeric(-fnew$var))
}

