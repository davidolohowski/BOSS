#' Two‐Stage BOSS Wrapper
#'
#' Runs the full two‐stage BOSS procedure:
#'   1) Bayesian sequential design via \code{BOSS_modal()} or \code{BOSS_mcmc()} (in development),
#'   2) Hessian update at the mode (only for \code{BOSS_modal()}),
#'   3) Essential‐support extraction and initial fill‐in,
#'   4) Approximate fill‐in to a target spacing,
#'   5) Final update of mode, Hessian, GP hyperparameters, and surrogate.
#'
#' @param func The objective function \code{f: [lower,upper]^D -> R}.
#' @param D Integer; input dimension.
#'
#' @param alpha Numeric in [0,1], significance level (default 0.05) for the essential support quantile.
#' The essential support is defined to capture around 1 - \code{alpha} of the posterior mass.
#' @param h Numeric; target fill-in distance (default 0.1).
#' @param method String; a character string specifying the BOSS method to use. Must be one of \code{moda} and \code{mcmc} (in development).
#'
#'
#' @param modal_opts List of options for \code{BOSS_modal()}. Any of
#'   \code{update_step}, \code{max_iter}, \code{quad}, \code{lower},
#'   \code{upper}, \code{nu}, \code{noise_var}, \code{modal_iter_check},
#'   \code{modal_check_warmup}, \code{modal_k.nn}, \code{modal_eps},
#'   \code{initial_design}, \code{delta}, \code{optim.n},
#'   \code{optim.max.iter}, \code{opt.lengthscale.grid}, \code{opt.grid}.
#' @param mcmc_opts (In development) List of options for \code{BOSS_mcmc()}. Any of
#' \code{update_step}, \code{max_iter}, \code{inner_iter}, \code{MCMC_size}, \code{quad},
#' \code{lower}, \code{upper}, \code{nu}, \code{noise_var}, \code{modal_k.nn},
#' \code{initial_design}, \code{delta}, \code{acc}, \code{explore_size},
#' \code{optim.n}, \code{optim.max.iter}, \code{Rhat_eps}, \code{UCB_prob}.
#' @param hess_opts List of options for \code{update_hessian()}. Any of
#'   \code{approach}, \code{GP.refine.eps}, \code{num.args}, \code{tol}.
#' @param essup_opts List of options for the essential‐support stage:
#'   \code{alpha} to construct the essential support and \code{n_samples}
#'   for \code{compute_fill_in()}.
#' @param fillin_opts List of options for the greedy \code{fill_in()} stage:
#'   \code{max_add}, \code{n_sample_max}.
#' @param verbose Integer 0–3; for progress.
#'
#' @return A fully‐updated S3 \code{boss} object with subclass being \code{boss_modal} or \code{boss_mcmc} depending on the argument of \code{method}.
#' @export
boss <- function(func,
                 D,
                 alpha = 0.05,
                 h = 0.1,
                 method = c('modal', 'mcmc'),
                 modal_opts  = list(),
                 mcmc_opts = list(),
                 hess_opts   = list(),
                 essup_opts  = list(),
                 fillin_opts = list(),
                 verbose     = 3) {

  method <- match.arg(method)

  if(!method %in% c('modal', 'mcmc')){
    stop("method must be one of 'modal' or 'mcmc'!")
  }

  if(method == 'modal' & length(mcmc_opts) > 0){
    warning(cat("BOSS modal method is supplied with arguments intended for BOSS mcmc. Make sure that the correct method was chosen."))
  }

  if(method == 'mcmc' & (length(modal_opts) > 0 | length(hess_opts) > 0)){
    warning(cat("BOSS mcmc method is supplied with arguments intended for BOSS modal. Make sure that the correct method was chosen."))
  }

  ## 1) Merge user options onto defaults
  modal_opts   <- setup_modal_opts(D, modal_opts)
  mcmc_opts    <- setup_mcmc_opts(D, alpha, mcmc_opts)
  hess_opts    <- setup_hess_opts(hess_opts)
  essup_opts   <- setup_essup_opts(alpha, essup_opts)
  fillin_opts  <- setup_fillin_opts(h, verbose, fillin_opts)

  if(method == 'modal'){
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
  }
  else{
    ## 2) Step 1: run BOSS_mcmc()
    br <- do.call(BOSS_mcmc,
                  c(list(func = func, D = D, verbose = verbose),
                    mcmc_opts))

    ## 3) Step 2: extract essential design points
    br <- construct_essential_designs(br)

    extra_essup <- essup_opts[setdiff(names(essup_opts),
                                      c("alpha"))]

    ## 4) compute_fill_in()
    br <- do.call(compute_fill_in,
                  c(list(boss_result = br), extra_essup))
  }
  ## 5) Step 4: greedy fill-in to target h
  if (verbose >= 1) {
    cat("Stage 2: Fill-in to target spacing h = ", fillin_opts$h, "\n")
    start_time <- Sys.time()
  }
  br <- do.call(fill_in,
                c(list(boss_result = br), fillin_opts))

  ## 6) Step 5: final update (mode, Hessian, GP hyperparams, surrogate)
  br <- update_boss(br, hess_opts = hess_opts)
  if (verbose >= 1) {
    end_time <- Sys.time()
    cat("Final update completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.\n")
  }

  return(br)
}





