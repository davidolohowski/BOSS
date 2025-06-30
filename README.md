
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BOSS

<!-- badges: start -->

<!-- badges: end -->

The goal of BOSS is to perform scalable approximate Bayesian inference
for complex hierarchical models, with a (log) marginal likelihood that
is intractable and/or expensive to compute.

## Installation

You can install the development version of BOSS from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("davidolohowski/BOSS")
```

## Example

This is a basic example which shows you how to solve a simple problem in
1-D.

Assume the true (log) marginal likelihood is given by the following
sinusoidal function:

``` r
library(BOSS)
f1d <- function(x) sin(0.5*pi*x)
```

We can use the `boss` function to approximate the (log) marginal
likelihood of this function, specifying the objective function, the
dimension of the problem, and the bounds for the BO optimization. The
`alpha` parameter controls the level of the essential support, while `h`
controls the target fill-in distance of the design points.

``` r
# A automatic implementation
br2 <- boss(f1d, D = 1, 
            modal_opts = list(lower = 0, upper = 2),
            alpha = 0.05, h = 0.1, verbose = 1)
#> Stage 1: Bayesian Optimization via Mode-Seeking Surrogate (BOSS) started.
#> Stage 1: BOSS finished.
#> Total time taken: 1.82 seconds.
#> Start updating Hessian at the mode...
#> Hessian updated in  0  seconds.
#> Stage 2: Fill-in to target spacing h =  0.1 
#> fill_in completed: added 12 point(s); final fill_in = 0.0651632 
#> Final update completed in  1.13  seconds.
plot(br2)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

For more advanced applications, the user can also manually perform each
step of the BOSS algorithm, as shown below:

``` r
# A manual implementation
br <- BOSS:::BOSS_modal(func = f1d, D = 1,
                        lower = 0, upper = 2,
                        update_step = 5,
                        max_iter = 20,
                        verbose = 1)
#> Stage 1: Bayesian Optimization via Mode-Seeking Surrogate (BOSS) started.
#> Stage 1: BOSS finished.
#> Total time taken: 0.93 seconds.
br <- update_hessian(br)
br <- compute_essential_support(br, alpha=0.05)
br <- construct_essential_designs(br)
br <- compute_fill_in(br)
br <- fill_in(br, h = 0.1, verbose = 1, max_add = 100)
#> fill_in completed: added 10 point(s); final fill_in = 0.097288
br_update <- update_boss(br)
plot(br_update)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />
