# This function computes the pairwise squared distances for matrices
compute_sq_dist <- function(X, Y) {
  if(is.vector(X)) X <- matrix(X, ncol = 1)
  if(is.vector(Y)) Y <- matrix(Y, ncol = 1)
  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  XY <- tcrossprod(X, Y)
  sq_dist <- matrix(XX, ncol=nrow(Y), nrow=nrow(X)) +
    t(matrix(YY, ncol=nrow(X), nrow=nrow(Y))) -
    2 * XY
  return(sq_dist)
}


# General covariance function generator: unified SE/Matérn via nu
cov_generator_nd <- function(length_scale = 1, signal_var = 1, nu = Inf) {

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
    x_s       <- sweep(x,        2, ls, "/")
    x_prime_s <- sweep(x_prime,  2, ls, "/")

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
      nu <- m + 0.5

      dist_s <- sqrt(sq_dist_s)
      sqrt_2nu <- sqrt(2 * nu)
      z <- sqrt_2nu * dist_s
      exp_part <- exp(-z)

      k_vec <- 0:m
      coeff_vec <- exp(lgamma(m + 1) - lgamma(2*m + 1) + lgamma(m + k_vec + 1) - (lgamma(k_vec + 1) + lgamma(m - k_vec + 1)))

      # Outer powers: z^k for each k
      z_vec <- as.vector(2*z)  # z_vec length n^2 if z is n x n

      # Vectorized: build z^k for each k (matrix n^2 x (m+1))
      z_powers <- outer(z_vec, m - k_vec, '^')  # fast

      # Sum across columns to get the final polynomial for each element
      poly_part_vec <- z_powers %*% coeff_vec

      # Reshape back to matrix
      poly_part <- matrix(poly_part_vec, nrow = nrow(z), ncol = ncol(z))

      cov_matrix <- signal_var * poly_part * exp_part
    }

    return(t(cov_matrix))
  }

  return(cov_fun)
}

#
# # Derivative generator of covariance matrix
# cov_deriv_generator_nd <- function(length_scale = 1, signal_var = 1, nu = Inf) {
#   cov_deriv_fun <- function(x, x_prime) {
#     # x_obs: n x D, x_star: 1 x D
#     if (is.vector(x))  x  <- matrix(x,  ncol = 1)
#     if (is.vector(x_prime)) x_prime <- matrix(x_prime, ncol = 1)
#     if (ncol(x) != ncol(x_prime))
#       stop("x and x_prime must have same number of columns")
#
#     n <- nrow(x)
#     D <- ncol(x)
#
#     # length‐scale vector
#     if (length(length_scale)==1) {
#       ls <- rep(length_scale, D)
#     } else {
#       ls <- length_scale
#     }
#     if (length(ls) != D)
#       stop("length_scale must be length 1 or match ncol(x)")
#
#     # scaled inputs
#     x_s      <- sweep(x,  2, ls, "/")    # n x D
#     x_prime_s <- sweep(x_prime, 2, ls, "/")    # 1 x D
#
#     # compute diffs and sq_dist
#     diffs   <- sweep(x_s, 2, as.numeric(x_prime_s), "-")  # n x D
#     sq_dist_s <- compute_sq_dist(x_s, x_prime_s)          # length n
#     sq_dist_s[sq_dist_s < 0] <- 0
#
#     # pre‐allocate derivative matrix
#     deriv <- matrix(0, n, D)
#
#     if (is.infinite(nu)) {
#       # --- SE kernel derivative ---
#       # base covariance: k = σ² exp(-½ sq_dist)
#       base_cov <- signal_var * exp(-sq_dist_s/2)           # length n
#       # ∂k/∂x'_j =  (x_j - x'_j)/ℓ² * k = diffs[,j]/ls[j] * base_cov
#       for (j in 1:D) {
#         deriv[,j] <- base_cov * ( diffs[,j] / ls[j] )
#       }
#
#     } else {
#       # --- half‐integer Matérn derivative ---
#       m   <- ceiling(nu) - 1
#       nu0 <- m + 0.5
#
#       dist_s <- sqrt(sq_dist_s)     # length n
#
#       sqrt2nu <- sqrt(2 * nu0)
#       z       <- sqrt2nu * dist_s         # length n
#       exp_z   <- exp(-z)                  # length n
#
#       # build P_m(z): polynomial of degree m
#       k_vec   <- 0:m
#       powers  <- m - k_vec
#       coeff_vec <- exp(lgamma(m + 1) - lgamma(2*m + 1) + lgamma(m + k_vec + 1) - (lgamma(k_vec + 1) + lgamma(m - k_vec + 1)))
#
#       # vectorized polynomial eval
#       # Outer powers: z^k for each k
#       z_vec <- as.vector(2*z)  # z_vec length n^2 if z is n x n
#
#       # Vectorized: build z^k for each k (matrix n^2 x (m+1))
#       z_powers <- outer(z_vec, powers, '^')
#       poly    <- z_powers %*% coeff_vec               # length n
#       poly <- matrix(poly, nrow = nrow(z), ncol = ncol(z))
#
#       # Build P_m'(z)
#       powers_prime <- pmax(powers - 1, 0)
#       coeff_prime <- coeff_vec * powers
#       z_powers_prime <- outer(z_vec, powers_prime, '^')
#       poly_prime <- z_powers_prime %*% coeff_prime
#       poly_prime <- matrix(poly_prime, nrow = nrow(z), ncol = ncol(z))
#
#       # dk/dr
#       dk_dr <- signal_var * sqrt2nu * exp_z * (2 * poly_prime - poly)
#
#       # finally ∂k/∂x'_j
#       for (j in 1:D) {
#         deriv[,j] <- -(diffs[,j] / dist_s / ls[j]) * dk_dr
#       }
#     }
#
#     return(deriv)  # n x D
#   }
#
#   return(cov_deriv_fun)
# }
#
#
# # Hessian generator of covariance matrix
# cov_hessian_generator_nd <- function(length_scale = 1, signal_var = 1, nu = Inf) {
#   cov_hessian_fun <- function(x, x_prime) {
#     # x: n x D, x_prime: 1 x D
#     if (is.vector(x)) x <- matrix(x, ncol = 1)
#     if (is.vector(x_prime)) x_prime <- matrix(x_prime, ncol = 1)
#     if (ncol(x) != ncol(x_prime))
#       stop("x and x_prime must have same number of columns")
#
#     n <- nrow(x)
#     D <- ncol(x)
#
#     # length‐scale vector
#     if (length(length_scale) == 1) {
#       ls <- rep(length_scale, D)
#     } else {
#       ls <- length_scale
#     }
#     if (length(ls) != D)
#       stop("length_scale must be length 1 or match ncol(x)")
#
#     # scaled inputs
#     x_s       <- sweep(x,       2, ls, "/")
#     x_prime_s <- sweep(x_prime, 2, ls, "/")
#
#     # compute diffs and sq_dist
#     diffs     <- sweep(x_s, 2, as.numeric(x_prime_s), "-")  # n x D
#     sq_dist_s <- compute_sq_dist(x_s, x_prime_s)            # length n
#     sq_dist_s[sq_dist_s < 0] <- 0
#
#     if (is.infinite(nu)) {
#       # --- SE kernel Hessian ---
#       base_cov <- signal_var * exp(-sq_dist_s/2)  # length n
#
#       H <- matrix(0, D, D)
#       for (i in 1:D) {
#         for (j in 1:D) {
#           term1 <- (diffs[,i] * diffs[,j]) / (ls[i]^2 * ls[j]^2)
#           delta_ij <- ifelse(i == j, 1, 0)
#           term2 <- delta_ij / (ls[i]^2)
#           H[i,j] <- sum(base_cov * (term1 - term2))
#         }
#       }
#
#     } else {
#       # --- half‐integer Matérn Hessian ---
#       m   <- ceiling(nu) - 1
#       nu0 <- m + 0.5
#
#       dist_s <- sqrt(sq_dist_s)   # length n
#       sqrt2nu <- sqrt(2 * nu0)
#       z       <- sqrt2nu * dist_s
#       exp_z   <- exp(-z)
#
#       # build P_m(z)
#       k_vec   <- 0:m
#       powers  <- m - k_vec
#       coeff_vec <- exp(lgamma(m + 1) - lgamma(2*m + 1) + lgamma(m + k_vec + 1) - (lgamma(k_vec + 1) + lgamma(m - k_vec + 1)))
#
#       z_vec <- as.vector(2*z)
#       z_powers <- outer(z_vec, powers, '^')
#       poly    <- z_powers %*% coeff_vec
#       poly <- matrix(poly, nrow = nrow(z), ncol = ncol(z))
#
#       # P_m'(z)
#       powers_prime <- pmax(powers - 1, 0)
#       coeff_prime <- coeff_vec * powers
#       z_powers_prime <- outer(z_vec, powers_prime, '^')
#       poly_prime <- z_powers_prime %*% coeff_prime
#       poly_prime <- matrix(poly_prime, nrow = nrow(z), ncol = ncol(z))
#
#       # P_m''(z)
#       powers_double_prime <- pmax(powers - 2, 0)
#       coeff_double_prime <- coeff_vec * powers * pmax(powers-1,0)
#       z_powers_double_prime <- outer(z_vec, powers_double_prime, '^')
#       poly_double_prime <- z_powers_double_prime %*% coeff_double_prime
#       poly_double_prime <- matrix(poly_double_prime, nrow = nrow(z), ncol = ncol(z))
#
#       dk_dr <- signal_var * sqrt2nu * exp_z * (2 * poly_prime - poly)
#       d2k_dr2 <- signal_var * (2 * nu0) * exp_z * (4 * poly_double_prime - 4 * poly_prime + poly)
#
#       H <- matrix(0, D, D)
#       for (i in 1:D) {
#         for (j in 1:D) {
#           delta_ij <- ifelse(i == j, 1, 0)
#           diffs_i <- diffs[,i]
#           diffs_j <- diffs[,j]
#
#           term1 <- (diffs_i * diffs_j) / (dist_s^2 * ls[i]^2 * ls[j]^2) * d2k_dr2
#           term2 <- (delta_ij / ls[i]^2 - (diffs_i * diffs_j) / (dist_s^2 * ls[i]^2 * ls[j]^2)) * dk_dr / dist_s
#
#           hess_contrib <- term1 + term2
#           hess_contrib[is.nan(hess_contrib)] <- 0  # handle 0/0 safely
#
#           # …then overwrite with the analytical Matérn limit at r=0
#           hess_contrib[dist_s == 0] <- if (i == j) - nu0/(nu0-1) * signal_var / (ls[i]^2) else 0
#
#           H[i,j] <- sum(hess_contrib)
#         }
#       }
#     }
#
#     return(H)  # D x D matrix
#   }
#
#   return(cov_hessian_fun)
# }

predict_gp <- function(data,
                       x_pred,
                       noise_var   = 1e-6,
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
    # — your original covariance + Cholesky —
    K_obs_obs  <- choice_cov(x_obs,  x_obs)
    K_obs_pred <- choice_cov(x_obs,  x_pred)
    K_pred_pred<- choice_cov(x_pred, x_pred)

    K_obs_obs <- K_obs_obs + noise_var * diag(N)
    L         <- chol(K_obs_obs)

    # solve for mean of f
    Ly        <- forwardsolve(t(L), y_obs)
    cond_mean <- K_obs_pred %*% backsolve(L, Ly)

    # solve for var of f
    LK        <- forwardsolve(t(L), t(K_obs_pred))
    cond_var  <- K_pred_pred - crossprod(LK)

    # — original single-GP fallback —
    sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)
    return(list(
      x    = x_pred,
      mean = as.vector(cond_mean),
      var  = cond_var,
      sim  = sim
    ))
  }

  else {
    # Compute covariance matrices
    K_obs_obs <- choice_cov(x_obs, x_obs)
    K_obs_pred <- choice_cov(x_obs, x_pred)
    K_pred_pred <- choice_cov(x_pred, x_pred)

    # Add noise to the diagonal and compute Cholesky decomposition
    K_obs_obs <- K_obs_obs + noise_var * diag(N)
    L <- chol(K_obs_obs)  # Upper triangular factor

    # Quadratic mean case
    X_covariate <- cbind(1, x_obs,
                         t(apply(x_obs, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))
    X_star <- cbind(1, x_pred,
                    t(apply(x_pred, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))

    # Solve using Cholesky decomposition
    LXs <- forwardsolve(t(L), X_covariate)
    #LX <- backsolve(L, LXs)
    A <- crossprod(LXs)
    Ly <- forwardsolve(t(L), y_obs)
    Ly <- backsolve(L, Ly)
    B <- crossprod(X_covariate, Ly)
    beta <- solve(A, B)

    # Compute residual solution
    residual <- y_obs - X_covariate %*% beta
    Lres <- forwardsolve(t(L), residual)
    Lres <- backsolve(L, Lres)
    cond_mean <- X_star %*% beta + K_obs_pred %*% Lres

    # Compute variance components
    LK_pred <- forwardsolve(t(L), t(K_obs_pred))
    var_part1 <- K_pred_pred - crossprod(LK_pred)

    delta <- X_star - t(crossprod(LXs, LK_pred))
    chol_A <- chol(A)
    Ldelta <- forwardsolve(t(chol_A), t(delta))
    var_part2 <- crossprod(Ldelta)
    cond_var <- var_part1 + var_part2

    # Generate simulation
    sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)

    return(list(x = x_pred, mean = as.vector(cond_mean), var = cond_var, sim = sim))
  }
}

# predict_gp <- function(data,
#                        x_pred,
#                        noise_var   = 1e-6,
#                        choice_cov,
#                        choice_cov_deriv   = NULL,
#                        choice_cov_hessian = NULL,
#                        quad = FALSE) {
#   # Extract x and y
#   x_obs <- data$x
#   y_obs <- data$y
#   if (is.vector(x_obs)) {
#     N <- length(x_obs); D <- 1
#   } else {
#     N <- nrow(x_obs);   D <- ncol(x_obs)
#   }
#
#   # — only in the non-quad + deriv+Hessian case do we do the multi-point joint —
#   if (!quad &&
#       !is.null(choice_cov_deriv) &&
#       !is.null(choice_cov_hessian) &&
#       nrow(x_pred) > 1) {
#
#     M <- nrow(x_pred)
#
#     # 1) build and Cholesky‐solve the training covariance
#     K_oo <- choice_cov(x_obs, x_obs)
#     K_oo <- K_oo + noise_var * diag(N)
#     L    <- chol(K_oo)  # upper‐triangular
#
#     # α = K_oo^{-1} y
#     v     <- forwardsolve(t(L), y_obs)
#     alpha <- backsolve(L, v)            # N×1
#
#     # 2) f* mean & var
#     K_op  <- choice_cov(x_obs, x_pred)   # N×M
#     W     <- forwardsolve(t(L), t(K_op))    # N×M
#     cond_mean <- as.vector(K_op %*% alpha)        # M×1
#     K_pp     <- choice_cov(x_pred, x_pred)           # M×M
#     cond_var <- K_pp - crossprod(W)                  # M×M
#
#     # 3) ∇f* mean & var
#     # 3a) build N×(M·D) matrix of cov(obs, deriv at each x*_i)
#     K_od_list <- lapply(1:M, function(i)
#       choice_cov_deriv(x_obs, x_pred[i,,drop=FALSE])  # N×D
#     )
#     K_od <- do.call(cbind, K_od_list)                 # N×(M·D)
#
#     # 3b) posterior mean of the stacked gradient vector
#     cond_deriv_mean <- as.vector(t(K_od) %*% alpha)    # (M·D)×1
#
#     K_grad_grad <- matrix(0, M*D, M*D)
#     for (i in seq_len(M)) {
#       for (j in seq_len(M)) {
#         H_ij <- choice_cov_hessian(
#           x_pred[i,,drop=FALSE],
#           x_pred[j,,drop=FALSE]
#         )  # D×D Hessian block
#         rows <- ((i-1)*D + 1):( i   *D)
#         cols <- ((j-1)*D + 1):( j   *D)
#         K_grad_grad[rows, cols] <- -H_ij
#       }
#     }
#
#     # 3c) posterior covariance of the stacked gradient
#     Z               <- forwardsolve(t(L), K_od)       # solves Lᵀ Z = K_od
#     Y               <- backsolve(      L, Z)          # solves L  Y = Z, so Y = K_oo⁻¹ K_od
#     cond_deriv_var  <- K_grad_grad - t(K_od) %*% Y     # (M·D)×(M·D)
#     # enforce exact symmetry
#     cond_deriv_var  <- (cond_deriv_var + t(cond_deriv_var)) / 2
#
#     # 4) prior cross‐cov between f* and ∇f*
#     #    build M×(M·D) block‐matrix
#     K_fp_list <- lapply(1:M, function(j)
#       choice_cov_deriv(x_pred, x_pred[j,,drop=FALSE])  # M×D
#     )
#     K_fp <- do.call(cbind, K_fp_list)                 # M×(M·D)
#
#     # —— FIX: reset any NaN (zero‐distance) entries to exactly zero
#     K_fp[is.nan(K_fp)] <- 0
#
#     cov_fgrad <- K_fp - t(W) %*% Z                  # M×(M·D)
#
#     # 6) assemble full joint covariance
#     top_left   <- cond_var                          #  M×M
#     top_right  <- cov_fgrad                         #  M×(M·D)
#     bot_left   <- t(cov_fgrad)                      # (M·D)×M
#     bot_right  <- cond_deriv_var                    # (M·D)×(M·D)
#     cov_joint  <- rbind(
#       cbind(top_left,  top_right),
#       cbind(bot_left,  bot_right)
#     )
#
#     # 7) joint simulate
#     mu_joint  <- c(cond_mean, cond_deriv_mean)
#     print(cov_joint)
#     sim_joint <- MASS::mvrnorm(1, mu_joint, cov_joint)
#
#     # split simulation
#     sim_f     <- sim_joint[     1:M]
#     sim_deriv <- matrix(sim_joint[-(1:M)], nrow=M, byrow=T)
#
#     return(list(
#       x            = x_pred,                     # M×D
#       mean         = cond_mean,                  # M×1
#       deriv_mean   = matrix(cond_deriv_mean, nrow=M, byrow = T),
#       var          = cond_var,                   # M×M
#       deriv_var    = cond_deriv_var,             # (M·D)×(M·D)
#       cov_joint    = cov_joint,                  # (M+M·D)×(M+M·D)
#       sim          = sim_f,                      # M×1
#       sim_deriv    = sim_deriv                   # M×D
#     ))
#   }
#
#   if (!quad) {
#     # — your original covariance + Cholesky —
#     K_obs_obs  <- choice_cov(x_obs,  x_obs)
#     K_obs_pred <- choice_cov(x_obs,  x_pred)
#     K_pred_pred<- choice_cov(x_pred, x_pred)
#
#     K_obs_obs <- K_obs_obs + noise_var * diag(N)
#     L         <- chol(K_obs_obs)
#
#     # solve for mean of f
#     Ly        <- forwardsolve(t(L), y_obs)
#     cond_mean <- K_obs_pred %*% backsolve(L, Ly)
#
#     # solve for var of f
#     LK        <- forwardsolve(t(L), t(K_obs_pred))
#     cond_var  <- K_pred_pred - crossprod(LK)
#
#     # — ADDED: if both deriv & hessian are supplied, do joint GP+∇GP —
#     if (!is.null(choice_cov_deriv) && !is.null(choice_cov_hessian)) {
#       # derivative cross-covariances
#       K_obs_pred_deriv  <- choice_cov_deriv(x_obs,  x_pred)    # N×D
#       K_pred_pred_deriv <- - choice_cov_hessian(x_pred, x_pred)  #   D×D
#
#       # posterior mean of ∇f
#       cond_deriv_mean <- t(K_obs_pred_deriv) %*% backsolve(L, Ly)  # D×1
#
#       # posterior var of ∇f
#       Z  <- forwardsolve(t(L), K_obs_pred_deriv)                   # N×D
#       Y  <- backsolve(L, Z)                                       # N×D = K⁻¹ K_obs_pred_deriv
#       cond_deriv_var <- K_pred_pred_deriv - t(K_obs_pred_deriv) %*% Y
#
#       # cross-covariance between f and ∇f
#       cov_fd <- -K_obs_pred %*% Y      # 1×D
#       cov_df <- t(cov_fd)                 # D×1
#
#       # assemble full joint covariance
#       top      <- cbind(cond_var,    cov_fd)       # 1×(1+D)
#       bottom   <- cbind(cov_df,    cond_deriv_var) # D×(1+D)
#       cov_joint<- rbind(top, bottom)               # (1+D)×(1+D)
#
#
#       # joint simulate
#       sim_joint <- MASS::mvrnorm(1,
#                                  c(as.vector(cond_mean),
#                                    as.vector(cond_deriv_mean)),
#                                  cov_joint)
#       sim        <- sim_joint[1]      # f
#       sim_deriv  <- sim_joint[-1]     # ∇f
#
#       return(list(
#         x            = x_pred,
#         mean         = as.vector(cond_mean),
#         deriv_mean   = as.vector(cond_deriv_mean),
#         var          = cond_var,
#         deriv_var    = cond_deriv_var,
#         cov_joint    = cov_joint,
#         sim          = sim,
#         sim_deriv    = sim_deriv
#       ))
#     }
#
#     # — original single-GP fallback —
#     sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)
#     return(list(
#       x    = x_pred,
#       mean = as.vector(cond_mean),
#       var  = cond_var,
#       sim  = sim
#     ))
#   }
#
#   else {
#     # Compute covariance matrices
#     K_obs_obs <- choice_cov(x_obs, x_obs)
#     K_obs_pred <- choice_cov(x_obs, x_pred)
#     K_pred_pred <- choice_cov(x_pred, x_pred)
#
#     # Add noise to the diagonal and compute Cholesky decomposition
#     K_obs_obs <- K_obs_obs + noise_var * diag(N)
#     L <- chol(K_obs_obs)  # Upper triangular factor
#
#     # Quadratic mean case
#     X_covariate <- cbind(1, x_obs,
#                          t(apply(x_obs, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))
#     X_star <- cbind(1, x_pred,
#                     t(apply(x_pred, 1, function(x) (x %o% x)[upper.tri(diag(D), TRUE)])))
#
#     # Solve using Cholesky decomposition
#     LXs <- forwardsolve(t(L), X_covariate)
#     #LX <- backsolve(L, LXs)
#     A <- crossprod(LXs)
#     Ly <- forwardsolve(t(L), y_obs)
#     Ly <- backsolve(L, Ly)
#     B <- crossprod(X_covariate, Ly)
#     beta <- solve(A, B)
#
#     # Compute residual solution
#     residual <- y_obs - X_covariate %*% beta
#     Lres <- forwardsolve(t(L), residual)
#     Lres <- backsolve(L, Lres)
#     cond_mean <- X_star %*% beta + K_obs_pred %*% Lres
#
#     # Compute variance components
#     LK_pred <- forwardsolve(t(L), t(K_obs_pred))
#     var_part1 <- K_pred_pred - crossprod(LK_pred)
#
#     delta <- X_star - t(crossprod(LXs, LK_pred))
#     chol_A <- chol(A)
#     Ldelta <- forwardsolve(t(chol_A), t(delta))
#     var_part2 <- crossprod(Ldelta)
#     cond_var <- var_part1 + var_part2
#
#     # Generate simulation
#     sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)
#
#     return(list(x = x_pred, mean = as.vector(cond_mean), var = cond_var, sim = sim))
#   }
# }


compute_like <- function(length_scale, x, y, signal_var, noise_var, D, prior_l_mean, prior_l_sd, quad = FALSE, nu = Inf){
  if(!quad){
    choice_cov <- cov_generator_nd(length_scale = length_scale, signal_var = signal_var, nu = nu)
    C <- choice_cov(x = x, x_prime = x)
    A <- C + noise_var * diag(nrow(x))
    L <- chol(A)  # A = L %*% t(L)

    # Compute the log-determinant of A and thus of Q = A^{-1}
    # (sum(log(diag(L))) = ½ log det(A))
    log_det_A <- sum(log(diag(L)))
    log_det_Q <- -log_det_A

    # Invert A via its Cholesky factor
    Q <- chol2inv(L)

    like <- as.numeric((t(y) %*% Q %*% y) / 2
                       - log_det_Q
                       - sum(dlnorm(length_scale, prior_l_mean + log(D)/2, prior_l_sd, log = TRUE)))
    if(is.nan(like) | is.na(like))      return(1e20)
    if(like == -Inf)                    return(-1e20)
    if(like == Inf)                     return(1e20)
    return(like)

  } else {
    ## --- set up design and covariance ---
    X_covariate <- matrix(cbind(rep(1, nrow(x)), x, t(apply(x, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))), nrow = nrow(x))

    choice_cov <- cov_generator_nd(length_scale = length_scale, signal_var = signal_var, nu = nu)
    C <- choice_cov(x = x, x_prime = x) + noise_var * diag(nrow(x))

    ## --- Cholesky-based inversion of C ---
    Lc <- chol(C)                  # C = t(Lc) %*% Lc  (R’s chol gives upper-triangular)
    Q  <- chol2inv(Lc)             # Q = C^{-1}

    ## --- build and invert the marginal precision S = X^T Q X via its Cholesky ---
    S_mat <- crossprod(X_covariate, Q %*% X_covariate)  # t(X_cov) %*% Q %*% X_cov
    Ls    <- chol(S_mat)            # S_mat = t(Ls) %*% Ls
    KK    <- chol2inv(Ls)           # KK = S_mat^{-1}

    ## --- log-determinants from the Cholesky factors ---
    #   sum(log(diag(Lc))) = ½ log det(C), so
    log_det_Q  <- -sum(log(diag(Lc)))  # = ½ log det(Q)
    #   sum(log(diag(Ls))) = ½ log det(S_mat), so
    log_det_KK <- -sum(log(diag(Ls)))  # = ½ log det(KK)

    ## --- compute the two quadratic forms ---
    q1 <- as.numeric((t(y) %*% Q %*% y) / 2)
    q2 <- as.numeric((t(y) %*% Q %*% X_covariate %*% KK %*%
                        t(X_covariate) %*% Q %*% y) / 2)

    ## --- assemble the log-likelihood (up to constants) ---
    like <- q1 - q2 - log_det_Q + log_det_KK -
      sum(dlnorm(length_scale, prior_l_mean + log(D)/2, prior_l_sd, log = TRUE))

    if(is.nan(like) | is.na(like))      return(1e20)
    if(like == -Inf)                    return(-1e20)
    if(like == Inf)                     return(1e20)
    return(like)
  }
}


UCB <- function(x, data, cov, nv, D, d, quad = FALSE){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv, quad = quad)

  # Compute the UCB acquisition function
  beta <- 2*log((D^2)*(pi^2)/(6*d))
  return(as.numeric(-fnew$mean - sqrt(beta*fnew$var)))
}

EXPLORE <- function(x, data, cov, nv, quad = FALSE){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv, quad = quad)

  # Compute the UCB acquisition function
  return(as.numeric(-fnew$var))
}


adaptive_MCMC_sampler <- function(N_sample = 15000, D, BO_surrogate, res_opt, prior){
  MCMC_sample <- matrix(0, nrow = N_sample, ncol = D)
  MCMC_sample[1,] <- res_opt$result$x[which.max(res_opt$result$y), ]
  for(i in 1:(N_sample - 1)) {
    x_old <- MCMC_sample[i, ]
    # Use a different proposal variance in early iterations
    if (i < 1000) {
      x_new <- x_old + rnorm(D, 0, 0.1)
    } else {
      # Adapt proposal using the empirical covariance from the chain
      prop_cov <- cov(MCMC_sample[1:i, , drop = FALSE]) + 1e-6 * diag(D)
      x_new <- MASS::mvrnorm(1, x_old, 2.38^2 / D * prop_cov)
    }
    # Define a local surrogate function for speed; use vectorized call inside optim-style evaluation
    f_val_new <- as.numeric(BO_surrogate(matrix(x_new, ncol = D)))
    f_val_old <- as.numeric(BO_surrogate(matrix(x_old, ncol = D)))
    # Proposal densities are uniform on [0,1] so add log-dunif values
    ratio <- f_val_new - f_val_old + sum(dunif(x_new, 0, 1, log = TRUE)) - sum(dunif(x_old, 0, 1, log = TRUE))
    if (log(runif(1)) < ratio && !is.na(ratio)) {
      MCMC_sample[i + 1, ] <- x_new
    } else {
      MCMC_sample[i + 1, ] <- x_old
    }
  }
  return(MCMC_sample)
}


sequential_BO_MCMC_modified <- function(func, lower, upper,
                                        control.list = list(max_iter = 50, inner_iter = 10, initial_design = 10,
                                                            alpha = 0.05, fill_number = 20,
                                                            GP = list(update_number = 10, nu = Inf, quad = FALSE,
                                                                      prior_l_mean = log(1), prior_l_sd = 0.5),
                                                            UCB = list(multistart_number = 1, delta = 0.01, k.nn = 5,
                                                                       acc = 0.5),
                                                            convergence = list(eps = 0.01, prob = 0.1))) {

  if(length(lower) != length(upper)){
    stop('Dimensions of lower and upper bound of the search domain does not match!')
  }

  D <- length(lower)

  # Store the original bounds (used for normalization in the GP)
  original_lower <- lower
  original_upper <- upper

  # Algorithm control
  max_iter <- control.list$max_iter # maximum iterations
  number_eval <- control.list$inner_iter # how many iterations to check for convergence
  initial_design <- control.list$initial_design # how many initial points to evaluate
  alpha <- control.list$alpha # What is quantile used for domain update
  fill_number <- control.list$fill_number # how many fill number to add at the end

  # GP control
  update_step <- control.list$GP$update_number # GP hyper-parameters update iteration
  nu <- control.list$GP$nu # GP covariance type (nu = Inf means SE; nu < Inf means matern)
  quad <- control.list$GP$quad # Do quadratic mean function or no
  prior_l_mean <- control.list$GP$prior_l_mean + log(D)/2 # prior for GP lengthscale
  prior_l_sd <- control.list$GP$prior_l_sd # prior for GP lengthscale

  # UCB control
  optim.n <- control.list$UCB$multistart_number # multistart number for UCB optimization
  delta <- control.list$UCB$delta # UCB delta parameter
  k.nn <- control.list$UCB$k.nn # number of nearest neighbor to check for modal performance
  acc <- control.list$UCB$acc # Level of accuracy required for modal convergence

  # Convergence control
  eps <- control.list$convergence$eps
  prob <- control.list$convergence$prob

  ## INITIAL DESIGN: Sample from the original region and normalize using original bounds
  xmat <- NULL          # Evaluations in original scale
  xmat_trans <- NULL    # Normalized evaluations (relative to original bounds)
  yvec <- c()

  cat("Initial random evaluation phase...\n")
  # initial <- matrix(runif(initial_design * D, original_lower, original_upper),
  #                   nrow = initial_design, ncol = D, byrow = TRUE)
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

  cat("  Updating GP hyperparameters...\n")

  if (!quad) {
    # no quadratic mean: just use the raw variance
    current_signal_var <- var(yvec)
  } else {
    # build the same quadratic design matrix that compute_like uses
    X_covariate <- matrix(cbind(rep(1, nrow(xmat_trans)),
                                xmat_trans,
                                t(apply(xmat_trans, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
                          nrow = nrow(xmat_trans))
    # fit the quadratic mean and get residuals
    fit_quad <- lm(yvec ~ X_covariate - 1)
    resids   <- resid(fit_quad)
    # signal variance is now the variance of those residuals
    current_signal_var <- var(resids)
  }

  opt <- optim(runif(D, 0.01, 0.9),
               function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
                                        signal_var = current_signal_var, noise_var = 1e-6,
                                        D, prior_l_mean, prior_l_sd, quad = quad, nu = nu),
               control = list(maxit = 100),
               lower = rep(0.01, D), upper = rep(0.9, D), method = "L-BFGS-B")
  current_length_scale <- opt$par
  lik <- opt$value
  cat("    GP log likelihood:", lik, "\n")
  cat("    New length_scale:", current_length_scale, "\n")
  cat("    New signal_var:", current_signal_var, "\n")

  # The "search region" used for acquisition and SPG injection.
  # Initially, this is the same as the original region.
  choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)

  cat("Initial MCMC run... \n")
  newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
  if (quad) {
    surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
      predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$sim
    }
  } else {
    surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
      predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$sim
    }
  }

  # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
  MCMC_sample <- adaptive_MCMC_sampler(N_sample = 15000, D = D, BO_surrogate = surrogate,
                                       res_opt = list(result = list(x = xmat_trans,
                                                                    x_original = xmat,
                                                                    y = yvec)))

  # Convert the normalized MCMC samples back to the original scale using original bounds.
  UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
  UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower

  explore_lower <- apply(MCMC_sample, 2, min) * (original_upper - original_lower) + original_lower
  explore_upper <- apply(MCMC_sample, 2, max) * (original_upper - original_lower) + original_lower

  cat("  New search region updated.\n")
  cat("  New lower bound:", UCB_lower, "\n")
  cat("  New upper bound:", UCB_upper, "\n")

  area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
  cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")

  # Update the search region used for candidate generation
  lower_current <- UCB_lower
  upper_current <- UCB_upper

  old_MCMC_sample <- MCMC_sample
  R_hat <- Inf

  UCB_counter <- 0

  ## MAIN BO LOOP (total BO iterations = max_iter)
  for (iter in 1:max_iter) {
    cat("BO Iteration:", iter, "\n")

    # Update GP hyperparameters every update_step iterations (using the GP data normalized by original bounds)
    if (iter %% update_step == 0) {
      cat("  Updating GP hyperparameters...\n")

      if (!quad) {
        # no quadratic mean: just use the raw variance
        current_signal_var <- var(yvec)
      } else {
        # build the same quadratic design matrix that compute_like uses
        X_covariate <- matrix(cbind(rep(1, nrow(xmat_trans)),
                                    xmat_trans,
                                    t(apply(xmat_trans, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
                              nrow = nrow(xmat_trans))
        # fit the quadratic mean and get residuals
        fit_quad <- lm(yvec ~ X_covariate - 1)
        resids   <- resid(fit_quad)
        # signal variance is now the variance of those residuals
        current_signal_var <- var(resids)
      }

      opt <- optim(runif(D, 0.01, 0.9),
                   function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
                                            signal_var = current_signal_var, noise_var = 1e-6,
                                            D, prior_l_mean, prior_l_sd, quad = quad, nu = nu),
                   control = list(maxit = 100),
                   lower = rep(0.01, D), upper = rep(0.9, D), method = "L-BFGS-B")
      current_length_scale <- opt$par
      lik <- opt$value
      cat("    GP log likelihood:", lik, "\n")
      cat("    New length_scale:", current_length_scale, "\n")
      cat("    New signal_var:", current_signal_var, "\n")
    }

    # Build the covariance function using the current hyperparameters
    choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)

    # Compute the adaptive probability for UCB based on the current design points.
    # Use the point with the best observed value (proxy for the mode).
    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
    if (quad) {
      surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
        predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
      }
    } else {
      surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
        predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
      }
    }


    # opt_check <- optim(init_point_trans,
    #                    function(x) surrogate(xvalue = matrix(x, nrow = 1, ncol = D)),
    #                    control = list(maxit = 100),
    #                    lower = (lower_current - original_lower) / (original_upper - original_lower),
    #                    upper = (upper_current - original_lower) / (original_upper - original_lower),
    #                    method = "L-BFGS-B", hessian = TRUE)
    mode_point <- xmat_trans[which.max(yvec),]
    mode_hess <- numDeriv::hessian(surrogate, matrix(mode_point, nrow = 1, ncol = D))
    mode_hess <- (mode_hess + t(mode_hess))/2
    c_const <- sqrt(sum(abs(eigen(mode_hess)$values)))/acc*log(3)
    # mode_point <- xmat_trans[which.max(yvec),]
    # Compute Euclidean distances from the best point to all other design points.
    distances <- sqrt(rowSums((xmat_trans - matrix(rep(mode_point, nrow(xmat_trans)),
                                                   nrow = nrow(xmat_trans), byrow = TRUE))^2))
    # Exclude the zero distance (self-distance)
    nonzero_distances <- distances[distances > 0]
    # If there is only one point, set delta_emp to a small value.
    # delta_emp <- if(length(nonzero_distances) == 0) 0 else median(nonzero_distances)
    d_min <- mean(nonzero_distances[order(nonzero_distances)[1:k.nn]])

    cat('  Current nearest neighbor distance around mode:', round(d_min*sqrt(sum((original_upper - original_lower)^2)), 6), '\n')

    p_value <- 2 / (1 + exp(-c_const * d_min)) - 1
    cat("  Adaptive UCB probability:", round(p_value, 4), "\n")

    UCB_flag <- runif(1) < p_value
    cat('  Do UCB?', UCB_flag, '\n')
    if(UCB_flag){
      lower_current <- UCB_lower
      upper_current <- UCB_upper

      # Maximize the acquisition function (here UCB) to select the next evaluation point.
      # For acquisition, we sample candidates from the current search region.
      # When normalizing for the GP, we always use the original bounds.
      # init_point <- matrix(runif(D*optim.n, lower_current, upper_current), nrow = optim.n, ncol = D, byrow = T)
      # init_point_trans <- t(apply(init_point, 1, function(x) (x - original_lower)/(original_upper - original_lower)))
      init_point <- lhs::randomLHS(optim.n, D)
      init_point <- t(apply(init_point, 1, function(x) x*(upper_current - lower_current) + lower_current))
      init_point_trans <- t(apply(init_point, 1, function(x) (x - original_lower)/(original_upper - original_lower)))

      cat("  Maximizing Acquisition Function...\n")
      optimizer <- optimx::multistart(init_point_trans,
                                function(x) UCB(matrix(x, nrow = 1, ncol = D),
                                                data = list(x = xmat_trans, x_original = xmat, y = yvec),
                                                cov = choice_cov, nv = 1e-6,
                                                D = UCB_counter + 1, d = delta, quad = quad),
                                control = list(maxit = 100),
                                lower = (lower_current - original_lower) / (original_upper - original_lower),
                                upper = (upper_current - original_lower) / (original_upper - original_lower),
                                method = "L-BFGS-B")

      next_point_trans <- unlist(unname(unique(optimizer[which.min(optimizer$value), c(1:D)])))
    }
    else{
      lower_current <- explore_lower
      upper_current <- explore_upper

      thinned <- unique(MCMC_sample[-c(1:1000),])
      s_thinned <- thinned[sample(1:nrow(thinned), 500),]

      cat("  Maximizing Acquisition Function...\n")
      var_thinned <- -apply(s_thinned, 1, function(x) EXPLORE(matrix(x, nrow = 1, ncol = D),
                                                              data = list(x = xmat_trans, x_original = xmat, y = yvec),
                                                              cov = choice_cov, nv = 1e-6, quad = quad))

      next_point_trans <- s_thinned[which.max(var_thinned),]

      # var_list <- parallel::mclapply(1:nrow(s_thinned), function(i) {
      #   # For each candidate row i, create a 1xD matrix and compute EXPLORE:
      #   EXPLORE(matrix(s_thinned[i, ], nrow = 1, ncol = D),
      #           data = list(x = xmat_trans, x_original = xmat, y = yvec),
      #           cov = choice_cov, nv = 1e-6, quad = quad)
      # }, mc.cores = num_cores)
      #
      # # Convert list to numeric vector and flip sign (as in your code)
      # var_thinned <- -unlist(var_list)
      # next_point_trans <- s_thinned[which.max(var_thinned), ]
    }

    UCB_counter <- UCB_counter + UCB_flag
    next_point <- next_point_trans * (original_upper - original_lower) + original_lower
    f_next <- func(next_point)

    # Append the new evaluation: for the GP we normalize using the original bounds.
    xmat_trans <- rbind(xmat_trans, (next_point - original_lower) / (original_upper - original_lower))
    xmat <- rbind(xmat, next_point)
    yvec <- c(yvec, f_next - rel)

    cat("  Next point:", paste(round(next_point, 4), collapse = ", "), "\n")
    cat("  Function value:", round(f_next, 4), "\n")

    ## MCMC UPDATE: Every number_eval BO iterations, update the search region via MCMC.
    if (iter %% number_eval == 0) {
      cat("  Running MCMC update at iteration", iter, "\n")
      newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
      if (quad) {
        surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
          predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$sim
        }
      } else {
        surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
          predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$sim
        }
      }

      # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
      MCMC_sample <- adaptive_MCMC_sampler(N_sample = 15000, D = D, BO_surrogate = surrogate,
                                           res_opt = list(result = list(x = xmat_trans,
                                                                        x_original = xmat,
                                                                        y = yvec)))

      # Convert the normalized MCMC samples back to the original scale using original bounds.
      UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
      UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower

      explore_lower <- apply(MCMC_sample, 2, min) * (original_upper - original_lower) + original_lower
      explore_upper <- apply(MCMC_sample, 2, max) * (original_upper - original_lower) + original_lower


      cat("  New search region updated.\n")
      cat("  New lower bound:", UCB_lower, "\n")
      cat("  New upper bound:", UCB_upper, "\n")

      area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
      cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")

      # Compute convergence diagnostic (R_hat)
      chains <- mcmc.list(as.mcmc(MCMC_sample[-c(1:1000), ]),
                          as.mcmc(old_MCMC_sample[-c(1:1000), ]))
      R_hat <- abs(max(coda::gelman.diag(chains)$psrf[, "Point est."]) - 1)
      cat("  R_hat of MCMC samples:", R_hat + 1, "\n")
      old_MCMC_sample <- MCMC_sample

      if (R_hat < eps & p_value < prob) {
        cat("Convergence achieved with R_hat =", round(R_hat, 5), "and UCB probability =", round(p_value, 4), "\n")
        break
      }
    }
  }

  if (R_hat >= eps | p_value >= prob) {
    cat("Maximum iterations reached and convergence criteria not met (R_hat =", round(R_hat, 5), "UCB probability =", round(p_value, 4), ")\n")
  }

  cat('  Inject sparse grid for filling... \n')


  samp_dat <- as.data.frame(MCMC_sample[-c(1:1000),])
  xmat_dat <- as.data.frame(xmat_trans)
  SPG <- as.matrix(dplyr::sample_n(dplyr::setdiff(samp_dat, xmat_dat), fill_number))

  SPG_original <- t(apply(SPG, 1, function(x) x*(original_upper - original_lower) + original_lower))

  # SPG <- t(apply(SPG_original, 1, function(x) (x - original_lower)/(original_upper - original_lower)))

  data_to_smooth <- list()
  unique_data <- unique(data.frame(x = xmat_trans, y = yvec + rel))
  data_to_smooth$x <- as.matrix(dplyr::select(unique_data, -y))
  data_to_smooth$x <- rbind(data_to_smooth$x, SPG)
  data_to_smooth$y <- unique_data$y

  cat('  Evaluate at sparse grid... \n')
  for(i in 1:nrow(SPG)){
    data_to_smooth$y <- c(data_to_smooth$y, func(SPG_original[i,]))
  }

  rel_y <- mean(data_to_smooth$y)
  data_to_smooth$y <- (data_to_smooth$y - rel_y)

  if (!quad) {
    # no quadratic mean: just use the raw variance
    current_signal_var <- var(data_to_smooth$y)
  } else {
    # build the same quadratic design matrix that compute_like uses
    X_covariate <- matrix(cbind(rep(1, nrow(data_to_smooth$x)),
                                data_to_smooth$x,
                                t(apply(data_to_smooth$x, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
                          nrow = nrow(data_to_smooth$x))
    # fit the quadratic mean and get residuals
    fit_quad <- lm(data_to_smooth$y ~ X_covariate - 1)
    resids   <- resid(fit_quad)
    # signal variance is now the variance of those residuals
    current_signal_var <- var(resids)
  }

  current_length_scale <- optim(runif(D, 0.01, 0.9), function(l) compute_like(length_scale = l, y = data_to_smooth$y, x = data_to_smooth$x,
                                                                              signal_var = current_signal_var, noise_var = 1e-6, D,
                                                                              prior_l_mean, prior_l_sd, quad = quad, nu = nu),
                                control = list(maxit = 100), lower = rep(0.01, D), upper = rep(0.9, D), method = 'L-BFGS-B')$par

  choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)

  # construct surrogate posterior
  if(quad){
    surrogate <- function(xvalue, dts = data_to_smooth, cov = choice_cov){
      predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
    }
  }
  else{
    surrogate <- function(xvalue, dts = data_to_smooth, cov = choice_cov){
      predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
    }
  }

  cat('  Final iteration of MCMC sampling from surrogate posterior... \n')
  final_MCMC_sample <- adaptive_MCMC_sampler(N_sample = 50000, D = D, BO_surrogate = surrogate,
                                             res_opt = list(result = list(x = xmat_trans,
                                                                          x_original = xmat,
                                                                          y = yvec)))

  MCMC_sample <- t(apply(final_MCMC_sample, 1, function(x) x*(original_upper - original_lower) + original_lower))

  res_opt <- list()
  res_opt$result$x <- data_to_smooth$x
  res_opt$result$x_original <- t(apply(data_to_smooth$x, 1, function(x) x*(original_upper - original_lower) + original_lower))
  res_opt$result$y <- data_to_smooth$y + rel_y
  res_opt$length_scale <- current_length_scale
  res_opt$signal_var <- current_signal_var

  return(list(BO_surrogate = res_opt, sample = MCMC_sample))
}


# MALA_sampler <- function(N_sample = 15000,
#                          D,
#                          BO_surrogate,
#                          prior,
#                          eps = 0.1) {
#
#   # storage
#   MCMC_sample <- matrix(0, nrow = N_sample, ncol = D)
#   # initialize at your BO optimum
#   MCMC_sample[1, ] <- rep(0.5, D)
#
#   # draw initial gradient at x0
#   out0      <- BO_surrogate(matrix(MCMC_sample[1, ], ncol = D))
#   grad_old  <- as.numeric(out0$sim_deriv)  # 1×D
#
#   for (i in seq_len(N_sample - 1)) {
#     x_old <- MCMC_sample[i, ]
#
#     # — MALA proposal drift + noise —
#     mu_prop <- x_old + (eps^2 / 2) * grad_old
#     x_new   <- mu_prop + eps * stats::rnorm(D)
#
#     # — single JOINT GP‐draw at both points —
#     pts <- rbind(x_old, x_new)  # 2×D
#     out <- BO_surrogate(pts)
#     f_pair    <- as.numeric(out$sim)        # length 2
#     grad_pair <- out$sim_deriv              # 2×D
#
#     f_old <- f_pair[1]
#     f_new <- f_pair[2]
#
#     # gradient at old & new under the same sample path
#     grad_old_est <- as.numeric(grad_pair[1, ])
#     grad_new     <- as.numeric(grad_pair[2, ])
#
#     # — log–target (surrogate + uniform prior) —
#     logp_old <- f_old + sum(dunif(x_old, 0, 1, log = TRUE))
#     logp_new <- f_new + sum(dunif(x_old, 0, 1, log = TRUE))
#
#     # — proposal densities —
#     logq_new_given_old <- sum(dnorm(x_new, mean = mu_prop, sd = eps, log = TRUE))
#     mu_rev            <- x_new + (eps^2 / 2) * grad_new
#     logq_old_given_new<- sum(dnorm(x_old, mean = mu_rev, sd = eps, log = TRUE))
#
#     # — MALA acceptance —
#     log_alpha <- logp_new - logp_old +
#       logq_old_given_new - logq_new_given_old
#
#     if (!is.na(log_alpha) && log(runif(1)) < log_alpha) {
#       MCMC_sample[i + 1, ] <- x_new
#       grad_old            <- grad_new       # for next iteration
#     } else {
#       MCMC_sample[i + 1, ] <- x_old
#       grad_old            <- grad_old_est   # update grad at x_old
#     }
#   }
#
#   return(MCMC_sample)
# }
#
#
# sequential_BO_MCMC_MALA <- function(func, lower, upper,
#                                         control.list = list(max_iter = 50, inner_iter = 10, initial_design = 10,
#                                                             alpha = 0.05, fill_number = 20,
#                                                             GP = list(update_number = 10, nu = Inf, quad = FALSE,
#                                                                       prior_l_mean = log(1), prior_l_sd = 0.5),
#                                                             UCB = list(multistart_number = 1, delta = 0.01, k.nn = 5,
#                                                                        acc = 0.5),
#                                                             convergence = list(eps = 0.01, prob = 0.1))) {
#
#   if(length(lower) != length(upper)){
#     stop('Dimensions of lower and upper bound of the search domain does not match!')
#   }
#
#   D <- length(lower)
#
#   # Store the original bounds (used for normalization in the GP)
#   original_lower <- lower
#   original_upper <- upper
#
#   # Algorithm control
#   max_iter <- control.list$max_iter # maximum iterations
#   number_eval <- control.list$inner_iter # how many iterations to check for convergence
#   initial_design <- control.list$initial_design # how many initial points to evaluate
#   alpha <- control.list$alpha # What is quantile used for domain update
#   fill_number <- control.list$fill_number # how many fill number to add at the end
#
#   # GP control
#   update_step <- control.list$GP$update_number # GP hyper-parameters update iteration
#   nu <- control.list$GP$nu # GP covariance type (nu = Inf means SE; nu < Inf means matern)
#   quad <- control.list$GP$quad # Do quadratic mean function or no
#   prior_l_mean <- control.list$GP$prior_l_mean + log(D)/2 # prior for GP lengthscale
#   prior_l_sd <- control.list$GP$prior_l_sd # prior for GP lengthscale
#
#   # UCB control
#   optim.n <- control.list$UCB$multistart_number # multistart number for UCB optimization
#   delta <- control.list$UCB$delta # UCB delta parameter
#   k.nn <- control.list$UCB$k.nn # number of nearest neighbor to check for modal performance
#   acc <- control.list$UCB$acc # Level of accuracy required for modal convergence
#
#   # Convergence control
#   eps <- control.list$convergence$eps
#   prob <- control.list$convergence$prob
#
#   ## INITIAL DESIGN: Sample from the original region and normalize using original bounds
#   xmat <- NULL          # Evaluations in original scale
#   xmat_trans <- NULL    # Normalized evaluations (relative to original bounds)
#   yvec <- c()
#
#   cat("Initial random evaluation phase...\n")
#   # initial <- matrix(runif(initial_design * D, original_lower, original_upper),
#   #                   nrow = initial_design, ncol = D, byrow = TRUE)
#   initial <- lhs::randomLHS(initial_design, D)
#   initial <- t(apply(initial, 1, function(x) x*(original_upper - original_lower) + original_lower))
#
#   for (i in 1:nrow(initial)) {
#     xmat_trans <- rbind(xmat_trans, (initial[i,] - original_lower) / (original_upper - original_lower))
#     xmat <- rbind(xmat, initial[i,])
#     yvec <- c(yvec, func(initial[i,]))
#   }
#
#   # Center the observations
#   rel <- mean(yvec)
#   yvec <- yvec - rel
#
#   cat("  Updating GP hyperparameters...\n")
#
#   if (!quad) {
#     # no quadratic mean: just use the raw variance
#     current_signal_var <- var(yvec)
#   } else {
#     # build the same quadratic design matrix that compute_like uses
#     X_covariate <- matrix(cbind(rep(1, nrow(xmat_trans)),
#                                 xmat_trans,
#                                 t(apply(xmat_trans, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
#                           nrow = nrow(xmat_trans))
#     # fit the quadratic mean and get residuals
#     fit_quad <- lm(yvec ~ X_covariate - 1)
#     resids   <- resid(fit_quad)
#     # signal variance is now the variance of those residuals
#     current_signal_var <- var(resids)
#   }
#
#   opt <- optim(runif(D, 0.01, 0.9),
#                function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
#                                         signal_var = current_signal_var, noise_var = 1e-6,
#                                         D, prior_l_mean, prior_l_sd, quad = quad, nu = nu),
#                control = list(maxit = 100),
#                lower = rep(0.01, D), upper = rep(0.9, D), method = "L-BFGS-B")
#   current_length_scale <- opt$par
#   lik <- opt$value
#   cat("    GP log likelihood:", lik, "\n")
#   cat("    New length_scale:", current_length_scale, "\n")
#   cat("    New signal_var:", current_signal_var, "\n")
#
#   # The "search region" used for acquisition and SPG injection.
#   # Initially, this is the same as the original region.
#   choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#   choice_cov_deriv <- cov_deriv_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#   choice_cov_hess <- cov_hessian_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#
#   cat("Initial MCMC run... \n")
#   newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
#   if (quad) {
#     surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#       predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#     }
#   } else {
#     surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#       predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#     }
#   }
#
#   # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
#   MCMC_sample <- MALA_sampler(N_sample = 15000, D = D, BO_surrogate = surrogate)
#
#   # Convert the normalized MCMC samples back to the original scale using original bounds.
#   UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
#   UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower
#
#   explore_lower <- apply(MCMC_sample, 2, min) * (original_upper - original_lower) + original_lower
#   explore_upper <- apply(MCMC_sample, 2, max) * (original_upper - original_lower) + original_lower
#
#   cat("  New search region updated.\n")
#   cat("  New lower bound:", UCB_lower, "\n")
#   cat("  New upper bound:", UCB_upper, "\n")
#
#   area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
#   cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")
#
#   # Update the search region used for candidate generation
#   lower_current <- UCB_lower
#   upper_current <- UCB_upper
#
#   old_MCMC_sample <- MCMC_sample
#   R_hat <- Inf
#
#   UCB_counter <- 0
#
#   ## MAIN BO LOOP (total BO iterations = max_iter)
#   for (iter in 1:max_iter) {
#     cat("BO Iteration:", iter, "\n")
#
#     # Update GP hyperparameters every update_step iterations (using the GP data normalized by original bounds)
#     if (iter %% update_step == 0) {
#       cat("  Updating GP hyperparameters...\n")
#
#       if (!quad) {
#         # no quadratic mean: just use the raw variance
#         current_signal_var <- var(yvec)
#       } else {
#         # build the same quadratic design matrix that compute_like uses
#         X_covariate <- matrix(cbind(rep(1, nrow(xmat_trans)),
#                                     xmat_trans,
#                                     t(apply(xmat_trans, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
#                               nrow = nrow(xmat_trans))
#         # fit the quadratic mean and get residuals
#         fit_quad <- lm(yvec ~ X_covariate - 1)
#         resids   <- resid(fit_quad)
#         # signal variance is now the variance of those residuals
#         current_signal_var <- var(resids)
#       }
#
#       opt <- optim(runif(D, 0.01, 0.9),
#                    function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
#                                             signal_var = current_signal_var, noise_var = 1e-6,
#                                             D, prior_l_mean, prior_l_sd, quad = quad, nu = nu),
#                    control = list(maxit = 100),
#                    lower = rep(0.01, D), upper = rep(0.9, D), method = "L-BFGS-B")
#       current_length_scale <- opt$par
#       lik <- opt$value
#       cat("    GP log likelihood:", lik, "\n")
#       cat("    New length_scale:", current_length_scale, "\n")
#       cat("    New signal_var:", current_signal_var, "\n")
#     }
#
#     # Build the covariance function using the current hyperparameters
#     choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#
#     # Compute the adaptive probability for UCB based on the current design points.
#     # Use the point with the best observed value (proxy for the mode).
#     newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
#     if (quad) {
#       surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
#         predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
#       }
#     } else {
#       surrogate <- function(xvalue, dts = newdata, cov = choice_cov) {
#         predict_gp(dts, x_pred = xvalue, choice_cov = cov, noise_var = 1e-6, quad = quad)$mean
#       }
#     }
#
#
#     # opt_check <- optim(init_point_trans,
#     #                    function(x) surrogate(xvalue = matrix(x, nrow = 1, ncol = D)),
#     #                    control = list(maxit = 100),
#     #                    lower = (lower_current - original_lower) / (original_upper - original_lower),
#     #                    upper = (upper_current - original_lower) / (original_upper - original_lower),
#     #                    method = "L-BFGS-B", hessian = TRUE)
#     mode_point <- xmat_trans[which.max(yvec),]
#     mode_hess <- numDeriv::hessian(surrogate, matrix(mode_point, nrow = 1, ncol = D))
#     mode_hess <- (mode_hess + t(mode_hess))/2
#     c_const <- sqrt(sum(abs(eigen(mode_hess)$values)))/acc*log(3)
#     # mode_point <- xmat_trans[which.max(yvec),]
#     # Compute Euclidean distances from the best point to all other design points.
#     distances <- sqrt(rowSums((xmat_trans - matrix(rep(mode_point, nrow(xmat_trans)),
#                                                    nrow = nrow(xmat_trans), byrow = TRUE))^2))
#     # Exclude the zero distance (self-distance)
#     nonzero_distances <- distances[distances > 0]
#     # If there is only one point, set delta_emp to a small value.
#     # delta_emp <- if(length(nonzero_distances) == 0) 0 else median(nonzero_distances)
#     d_min <- mean(nonzero_distances[order(nonzero_distances)[1:k.nn]])
#
#     cat('  Current nearest neighbor distance around mode:', round(d_min*sqrt(sum((original_upper - original_lower)^2)), 6), '\n')
#
#     p_value <- 2 / (1 + exp(-c_const * d_min)) - 1
#     cat("  Adaptive UCB probability:", round(p_value, 4), "\n")
#
#     UCB_flag <- runif(1) < p_value
#     cat('  Do UCB?', UCB_flag, '\n')
#     if(UCB_flag){
#       lower_current <- UCB_lower
#       upper_current <- UCB_upper
#
#       # Maximize the acquisition function (here UCB) to select the next evaluation point.
#       # For acquisition, we sample candidates from the current search region.
#       # When normalizing for the GP, we always use the original bounds.
#       # init_point <- matrix(runif(D*optim.n, lower_current, upper_current), nrow = optim.n, ncol = D, byrow = T)
#       # init_point_trans <- t(apply(init_point, 1, function(x) (x - original_lower)/(original_upper - original_lower)))
#       init_point <- lhs::randomLHS(optim.n, D)
#       init_point <- t(apply(init_point, 1, function(x) x*(upper_current - lower_current) + lower_current))
#       init_point_trans <- t(apply(init_point, 1, function(x) (x - original_lower)/(original_upper - original_lower)))
#
#       cat("  Maximizing Acquisition Function...\n")
#       optimizer <- optimx::multistart(init_point_trans,
#                                       function(x) UCB(matrix(x, nrow = 1, ncol = D),
#                                                       data = list(x = xmat_trans, x_original = xmat, y = yvec),
#                                                       cov = choice_cov, nv = 1e-6,
#                                                       D = UCB_counter + 1, d = delta, quad = quad),
#                                       control = list(maxit = 100),
#                                       lower = (lower_current - original_lower) / (original_upper - original_lower),
#                                       upper = (upper_current - original_lower) / (original_upper - original_lower),
#                                       method = "L-BFGS-B")
#
#       next_point_trans <- unlist(unname(unique(optimizer[which.min(optimizer$value), c(1:D)])))
#     }
#     else{
#       lower_current <- explore_lower
#       upper_current <- explore_upper
#
#       thinned <- unique(MCMC_sample[-c(1:1000),])
#       s_thinned <- thinned[sample(1:nrow(thinned), 500),]
#
#       cat("  Maximizing Acquisition Function...\n")
#       var_thinned <- -apply(s_thinned, 1, function(x) EXPLORE(matrix(x, nrow = 1, ncol = D),
#                                                               data = list(x = xmat_trans, x_original = xmat, y = yvec),
#                                                               cov = choice_cov, nv = 1e-6, quad = quad))
#
#       next_point_trans <- s_thinned[which.max(var_thinned),]
#
#       # var_list <- parallel::mclapply(1:nrow(s_thinned), function(i) {
#       #   # For each candidate row i, create a 1xD matrix and compute EXPLORE:
#       #   EXPLORE(matrix(s_thinned[i, ], nrow = 1, ncol = D),
#       #           data = list(x = xmat_trans, x_original = xmat, y = yvec),
#       #           cov = choice_cov, nv = 1e-6, quad = quad)
#       # }, mc.cores = num_cores)
#       #
#       # # Convert list to numeric vector and flip sign (as in your code)
#       # var_thinned <- -unlist(var_list)
#       # next_point_trans <- s_thinned[which.max(var_thinned), ]
#     }
#
#     UCB_counter <- UCB_counter + UCB_flag
#     next_point <- next_point_trans * (original_upper - original_lower) + original_lower
#     f_next <- func(next_point)
#
#     # Append the new evaluation: for the GP we normalize using the original bounds.
#     xmat_trans <- rbind(xmat_trans, (next_point - original_lower) / (original_upper - original_lower))
#     xmat <- rbind(xmat, next_point)
#     yvec <- c(yvec, f_next - rel)
#
#     cat("  Next point:", paste(round(next_point, 4), collapse = ", "), "\n")
#     cat("  Function value:", round(f_next, 4), "\n")
#
#     ## MCMC UPDATE: Every number_eval BO iterations, update the search region via MCMC.
#     if (iter %% number_eval == 0) {
#       cat("  Running MCMC update at iteration", iter, "\n")
#       newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
#       if (quad) {
#         surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#           predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#         }
#       } else {
#         surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#           predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#         }
#       }
#
#       # Run the adaptive MCMC sampler (note: it works in normalized [0,1] space)
#       MCMC_sample <- MALA_sampler(N_sample = 15000, D = D, BO_surrogate = surrogate)
#
#       # Convert the normalized MCMC samples back to the original scale using original bounds.
#       UCB_lower <- apply(MCMC_sample, 2, quantile, probs = alpha) * (original_upper - original_lower) + original_lower
#       UCB_upper <- apply(MCMC_sample, 2, quantile, probs = 1-alpha) * (original_upper - original_lower) + original_lower
#
#       explore_lower <- apply(MCMC_sample, 2, min) * (original_upper - original_lower) + original_lower
#       explore_upper <- apply(MCMC_sample, 2, max) * (original_upper - original_lower) + original_lower
#
#
#       cat("  New search region updated.\n")
#       cat("  New lower bound:", UCB_lower, "\n")
#       cat("  New upper bound:", UCB_upper, "\n")
#
#       area <- prod(UCB_upper - UCB_lower) / prod(original_upper - original_lower)
#       cat("  New search region contracted by", round((1 - area) * 100, 2), "%\n")
#
#       # Compute convergence diagnostic (R_hat)
#       chains <- mcmc.list(as.mcmc(MCMC_sample[-c(1:1000), ]),
#                           as.mcmc(old_MCMC_sample[-c(1:1000), ]))
#       R_hat <- abs(max(coda::gelman.diag(chains)$psrf[, "Point est."]) - 1)
#       cat("  R_hat of MCMC samples:", R_hat + 1, "\n")
#       old_MCMC_sample <- MCMC_sample
#
#       if (R_hat < eps & p_value < prob) {
#         cat("Convergence achieved with R_hat =", round(R_hat, 5), "and UCB probability =", round(p_value, 4), "\n")
#         break
#       }
#     }
#   }
#
#   if (R_hat >= eps | p_value >= prob) {
#     cat("Maximum iterations reached and convergence criteria not met (R_hat =", round(R_hat, 5), "UCB probability =", round(p_value, 4), ")\n")
#   }
#
#   cat('  Inject sparse grid for filling... \n')
#
#
#   samp_dat <- as.data.frame(MCMC_sample[-c(1:1000),])
#   xmat_dat <- as.data.frame(xmat_trans)
#   SPG <- as.matrix(dplyr::sample_n(dplyr::setdiff(samp_dat, xmat_dat), fill_number))
#
#   SPG_original <- t(apply(SPG, 1, function(x) x*(original_upper - original_lower) + original_lower))
#
#   # SPG <- t(apply(SPG_original, 1, function(x) (x - original_lower)/(original_upper - original_lower)))
#
#   data_to_smooth <- list()
#   unique_data <- unique(data.frame(x = xmat_trans, y = yvec + rel))
#   data_to_smooth$x <- as.matrix(dplyr::select(unique_data, -y))
#   data_to_smooth$x <- rbind(data_to_smooth$x, SPG)
#   data_to_smooth$y <- unique_data$y
#
#   cat('  Evaluate at sparse grid... \n')
#   for(i in 1:nrow(SPG)){
#     data_to_smooth$y <- c(data_to_smooth$y, func(SPG_original[i,]))
#   }
#
#   rel_y <- mean(data_to_smooth$y)
#   data_to_smooth$y <- (data_to_smooth$y - rel_y)
#
#   if (!quad) {
#     # no quadratic mean: just use the raw variance
#     current_signal_var <- var(data_to_smooth$y)
#   } else {
#     # build the same quadratic design matrix that compute_like uses
#     X_covariate <- matrix(cbind(rep(1, nrow(data_to_smooth$x)),
#                                 data_to_smooth$x,
#                                 t(apply(data_to_smooth$x, 1, function(y) (((y%*%t(y))[upper.tri(y%*%t(y), diag = T)]))))),
#                           nrow = nrow(data_to_smooth$x))
#     # fit the quadratic mean and get residuals
#     fit_quad <- lm(data_to_smooth$y ~ X_covariate - 1)
#     resids   <- resid(fit_quad)
#     # signal variance is now the variance of those residuals
#     current_signal_var <- var(resids)
#   }
#
#   current_length_scale <- optim(runif(D, 0.01, 0.9), function(l) compute_like(length_scale = l, y = data_to_smooth$y, x = data_to_smooth$x,
#                                                                               signal_var = current_signal_var, noise_var = 1e-6, D,
#                                                                               prior_l_mean, prior_l_sd, quad = quad, nu = nu),
#                                 control = list(maxit = 100), lower = rep(0.01, D), upper = rep(0.9, D), method = 'L-BFGS-B')$par
#
#   choice_cov <- cov_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#   choice_cov_deriv <- cov_deriv_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#   choice_cov_hess <- cov_hessian_generator_nd(length_scale = current_length_scale, signal_var = current_signal_var, nu = nu)
#
#   # construct surrogate posterior
#   if (quad) {
#     surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#       predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#     }
#   } else {
#     surrogate <- function(xvalue, dts = newdata, cov = choice_cov, deriv = choice_cov_deriv, hess = choice_cov_hess) {
#       predict_gp(dts, x_pred = xvalue, choice_cov = cov, choice_cov_deriv = deriv, choice_cov_hessian = hess, noise_var = 1e-6, quad = quad)
#     }
#   }
#
#   cat('  Final iteration of MCMC sampling from surrogate posterior... \n')
#   final_MCMC_sample <- MALA_sampler(N_sample = 50000, D = D, BO_surrogate = surrogate)
#
#   MCMC_sample <- t(apply(final_MCMC_sample, 1, function(x) x*(original_upper - original_lower) + original_lower))
#
#   res_opt <- list()
#   res_opt$result$x <- data_to_smooth$x
#   res_opt$result$x_original <- t(apply(data_to_smooth$x, 1, function(x) x*(original_upper - original_lower) + original_lower))
#   res_opt$result$y <- data_to_smooth$y + rel_y
#   res_opt$length_scale <- current_length_scale
#   res_opt$signal_var <- current_signal_var
#
#   return(list(BO_surrogate = res_opt, sample = MCMC_sample))
# }
#
#
#
#
#
