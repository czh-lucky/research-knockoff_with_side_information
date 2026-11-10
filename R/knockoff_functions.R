create.gaussian.fast <- function (X, mu, Sigma, method = c("asdp", "sdp", "equi"), diag_s = NULL) 
{
  method = match.arg(method)
  if ((nrow(Sigma) <= 500) && method == "asdp") {
    method = "sdp"
  }
  if (is.null(diag_s)) {
    diag_s = diag(switch(match.arg(method), 
                         equi = knockoff:::create.solve_equi(Sigma), 
                         sdp = knockoff:::create.solve_sdp(Sigma), 
                         asdp = create.solve_asdp.fast(Sigma)))
  }
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s, length(diag_s))
  }
  if (all(diag_s == 0)) {
    warning("The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power.")
    return(X)
  }

  SigmaInv_s = solve(Sigma, diag_s)

  mu_k = X - sweep(X, 2, mu, "-") %*% SigmaInv_s

  s <- diag(diag_s)
  Sigma_k <- -sweep(SigmaInv_s, 1, s, "*")
  diag(Sigma_k) <- diag(Sigma_k) + 2 * s

  X_k = mu_k + matrix(
    rnorm(ncol(X) * nrow(X)),
    nrow(X)
  ) %*% chol(Sigma_k)
}

create.second_order.fast <- function (
    X,
    method = c("asdp", "equi", "sdp"),
    shrink = F
) 
{
  method = match.arg(method)

  mu = colMeans(X)

  if (!shrink) {
    Sigma = cov(X)

    if (!knockoff:::is_posdef(Sigma)) {
      shrink = TRUE
    }
  }

  if (shrink) {
    if (!requireNamespace("corpcor", quietly = T)) 
      stop("corpcor is not installed", call. = F)

    Sigma = tryCatch({
      suppressWarnings(
        matrix(
          as.numeric(
            corpcor::cov.shrink(X, verbose = F)
          ),
          nrow = ncol(X)
        )
      )
    }, warning = function(w) {
    }, error = function(e) {
      stop(
        "SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0"
      )
    }, finally = {
    })
  }

  create.gaussian.fast(X, mu, Sigma, method = method)
}


create.solve_asdp.fast <- function(
    Sigma,
    max.size = 500,
    gaptol = 1e-6,
    maxit = 1000,
    verbose = FALSE
) {
    stopifnot(isSymmetric(Sigma))

    if (ncol(Sigma) <= max.size) {
        return(
            knockoff:::create.solve_sdp(
                Sigma,
                gaptol = gaptol,
                maxit = maxit,
                verbose = verbose
            )
        )
    }

    if (verbose) {
        cat(sprintf(
            "Dividing the problem into subproblems of size <= %s ... ",
            max.size
        ))
    }

    cluster_sol <- knockoff:::divide.sdp(
        Sigma,
        max.size = max.size
    )

    n.blocks <- max(cluster_sol$clusters)

    if (verbose)
        cat("done.\n")

    if (verbose) {
        cat(sprintf(
            "Solving %s smaller SDPs ...\n",
            n.blocks
        ))
    }

    s_asdp_list <- list()

    if (verbose) {
        pb <- utils::txtProgressBar(
            min = 0,
            max = n.blocks,
            style = 3
        )
    }

    for (k in 1:n.blocks) {
        s_asdp_list[[k]] <- knockoff:::create.solve_sdp(
            as.matrix(cluster_sol$subSigma[[k]]),
            gaptol = gaptol,
            maxit = maxit
        )

        if (verbose)
            utils::setTxtProgressBar(pb, k)
    }

    if (verbose)
        cat("\n")

    p <- dim(Sigma)[1]

    idx_count <- rep(1, n.blocks)
    s_asdp <- rep(0, p)

    for (j in 1:p) {
        cluster_j <- cluster_sol$clusters[j]

        s_asdp[j] <-
            s_asdp_list[[cluster_j]][
                idx_count[cluster_j]
            ]

        idx_count[cluster_j] <-
            idx_count[cluster_j] + 1
    }

    if (verbose) {
        cat(sprintf(
            "Combining the solutions of the %s smaller SDPs ... ",
            n.blocks
        ))
    }

    cat("\nUnscaled s_asdp summary:\n")
    print(summary(s_asdp))

    ## Gamma optimization

    if (all(s_asdp == 0)) {
        if (verbose)
            cat("done.\n")

        return(s_asdp)
    }

    ## chol() returns upper triangular R:
    ##
    ## Sigma = t(R) %*% R
    ##
    R <- chol(Sigma)

    ## Matrix-free multiplication by
    ##
    ## B = R^{-T} D R^{-1},
    ##
    ## where D = diag(s_asdp).
    B_mult <- function(x, args) {
        R <- args$R
        s <- args$s

        ## y = R^{-1} x
        y <- backsolve(R, x)

        ## z = D y
        z <- s * y

        ## w = R^{-T} z
        w <- backsolve(
            R,
            z,
            transpose = TRUE
        )

        as.numeric(w)
    }

    ## Compute lambda_max(B)
    eig_out <- RSpectra::eigs_sym(
        B_mult,
        k = 1,
        which = "LA",
        n = p,
        args = list(
            R = R,
            s = s_asdp
        ),
        opts = list(
            retvec = FALSE,
            maxitr = 1e5,
            tol = 1e-12
        )
    )

    lambda_max <- as.numeric(
        eig_out$values[1]
    )

    if (!is.finite(lambda_max)) {
        stop(
            "Failed to compute the largest eigenvalue in gamma optimization."
        )
    }

    if (lambda_max <= 0) {
        gamma_bound <- Inf
    } else {
        gamma_bound <- 2 / lambda_max
    }

    ## Project gamma back to the original 1000-point grid.
    ## Use the largest grid point strictly below gamma_bound.
    gamma_range <- seq(
        0,
        1,
        length.out = 1000
    )

    if (gamma_bound > 1) {

        gamma <- 1

    } else {

        valid_gamma <- gamma_range[
            gamma_range < gamma_bound
        ]

        if (length(valid_gamma) == 0) {
            gamma <- 0
        } else {
            gamma <- max(valid_gamma)
        }
    }

    if (verbose)
        cat("Verifying that the solution is correct ... ")
    
    s_asdp_scaled <- gamma * s_asdp

    valid <- knockoff:::is_posdef(
        2 * Sigma - diag(s_asdp_scaled, length(s_asdp_scaled))
    )

    if (!valid) {
        gamma_idx <- which.min(abs(gamma_range - gamma))
        if (gamma_idx > 1){
            gamma <- gamma_range[gamma_idx - 1]
            s_asdp_scaled <- gamma * s_asdp

            valid <- knockoff:::is_posdef(
                2 * Sigma - diag(s_asdp_scaled, length(s_asdp_scaled))
            )
        }
    }

    if (!valid){
        warning(
            paste(
                "In creation of approximate SDP knockoffs,",
                "procedure failed. Knockoffs will have no power."
            ),
            immediate. = TRUE
        )

        s_asdp_scaled <- 0 * s_asdp_scaled
    }

    cat("\nScaled s summary:\n")
    print(summary(s_asdp_scaled))

    cat(sprintf(
        "lambda_max(B) = %.12g\n",
        lambda_max
    ))

    cat(sprintf(
        "gamma_bound   = %.12g\n",
        gamma_bound
    ))

    cat(sprintf(
        "gamma         = %.12g\n",
        gamma
    ))

    cat(sprintf(
        "positive definite = %s\n", valid
    ))
    
    if (verbose)
        cat("done.\n")

    s_asdp_scaled
}