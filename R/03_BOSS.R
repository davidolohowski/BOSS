#' Two‐Stage BOSS Wrapper
#'
#' Runs the full two‐stage BOSS procedure:
#'   1) Bayesian sequential design via \code{BOSS_modal()},
#'   2) Hessian update at the mode,
#'   3) Essential‐support extraction and initial fill‐in,
#'   4) Approximate fill‐in to a target spacing,
#'   5) Final update of mode, Hessian, GP hyperparameters, and surrogate.
#'
#' @param func The objective function \code{f: [lower,upper]^D -> R}.
#' @param D Integer; input dimension.
#'
#' @param alpha Numeric in [0,1], significance level (default 0.05) for the chi-squared cutoff.
#' The essential support is defined to capture around 1 - \code{alpha} of the posterior mass.
#' @param h Numeric; target fill-in distance (default 0.1).
#'
#'
#' @param modal_opts List of options for \code{BOSS_modal()}. Any of
#'   \code{update_step}, \code{max_iter}, \code{quad}, \code{lower},
#'   \code{upper}, \code{nu}, \code{noise_var}, \code{modal_iter_check},
#'   \code{modal_check_warmup}, \code{modal_k.nn}, \code{modal_eps},
#'   \code{initial_design}, \code{delta}, \code{optim.n},
#'   \code{optim.max.iter}, \code{opt.lengthscale.grid}, \code{opt.grid}.
#' @param hess_opts List of options for \code{update_hessian()}. Any of
#'   \code{approach}, \code{local.poly.args}, \code{num.args}, \code{tol}.
#' @param essup_opts List of options for the essential‐support stage:
#'   \code{alpha} to construct the essential support and \code{n_samples}
#'   for \code{compute_fill_in()}.
#' @param fillin_opts List of options for the greedy \code{fill_in()} stage:
#'   \code{max_add}, \code{n_sample_max}.
#' @param verbose Integer 0–3; for progress.
#'
#' @return A fully‐updated S3 \code{boss} object.
#' @export
boss <- function(func,
                 D,
                 alpha = 0.05,
                 h = 0.1,
                 modal_opts  = list(),
                 hess_opts   = list(),
                 essup_opts  = list(),
                 fillin_opts = list(),
                 verbose     = 3) {

  ## 0) Default option lists
  default_modal <- list(
    update_step         = 5,
    max_iter            = 100,
    quad                = FALSE,
    lower               = rep(0, D),
    upper               = rep(1, D),
    nu                  = Inf,
    noise_var           = 1e-6,
    modal_iter_check    = 10,
    modal_check_warmup  = 20,
    modal_k.nn          = 5,
    modal_eps           = 0.1,
    initial_design      = 5,
    delta               = 0.01,
    optim.n             = 5,
    optim.max.iter      = 1000,
    opt.lengthscale.grid= NULL,
    opt.grid            = NULL
  )
  default_hess  <- list(
    approach      = 'local.poly',
    local.poly.args = list(eps = 0.1, bw = NULL, kernel = 'Epanech'),
    num.args = list(method = 'Richardson', method.args = list()),
    tol         = 1e-8
  )
  default_essup <- list(
    alpha             = alpha,
    n_samples         = 10000
  )
  default_fillin <- list(
    h         = h,
    max_add   = 100,
    n_sample_max = 10000,
    verbose   = verbose
  )

  ## 1) Merge user options onto defaults
  modal_opts   <- modifyList(default_modal,   modal_opts)
  hess_opts    <- modifyList(default_hess,    hess_opts)
  essup_opts   <- modifyList(default_essup,   essup_opts)
  fillin_opts  <- modifyList(default_fillin,  fillin_opts)

  ## 2) Step 1: run BOSS_modal()
  br <- do.call(BOSS_modal,
                c(list(func = func, D = D, verbose = verbose),
                  modal_opts))

  ## 3) Step 2: update Hessian at the mode
  if (verbose >= 1) {
    cat("Start updating Hessian at the mode...\n")
    start_time <- Sys.time()
  }
  br <- do.call(update_hessian,
                c(list(boss_result = br),
                  hess_opts))
  if (verbose >= 1) {
    end_time <- Sys.time()
    cat("Hessian updated in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.\n")
  }

  ## 4) Step 3: essential‐support + initial fill‐in
  br <- compute_essential_support(br, alpha = essup_opts$alpha)
  br <- construct_essential_designs(br)
  # pass the remaining essup arguments to compute_fill_in()
  extra_essup <- essup_opts[setdiff(names(essup_opts),
                                    c("alpha"))]
  br <- do.call(compute_fill_in,
                c(list(boss_result = br), extra_essup))

  ## 5) Step 4: greedy fill-in to target h
  if (verbose >= 1) {
    cat("Stage 2: Fill-in to target spacing h = ", fillin_opts$h, "\n")
    start_time <- Sys.time()
  }
  br <- do.call(fill_in,
                c(list(boss_result = br), fillin_opts))

  ## 6) Step 5: final update (mode, Hessian, GP hyperparams, surrogate)
  br <- update_boss(br)
  if (verbose >= 1) {
    end_time <- Sys.time()
    cat("Final update completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.\n")
  }

  return(br)
}





