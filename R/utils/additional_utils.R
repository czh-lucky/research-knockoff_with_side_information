filter_gam <- function(W, U, fdr = 0.1, offset = 1,
                       df = 5, reveal_prop = 0.5, mute = TRUE){

  ## Determine dimensionality
  p <- length(W)
  
  ## Initilization
  tau <- rep(0,p)
  rejs <- list()
  W_abs <- abs(W)
  W_sign <- as.numeric(W>0)
  revealed_sign <- rep(1,p)
  all_id <- 1:p
  tau.sel <- c()

  ## Reveal the features with |W|>=quantile(|W|,reveal_prop)
  revealed_id <- which(abs(W) <= quantile(abs(W), reveal_prop))

  ## Update the revealed information
  revealed_sign[revealed_id] <- W_sign[revealed_id]
  unrevealed_id <- all_id[-revealed_id]
  tau[revealed_id] <- W_abs[revealed_id] + 1

  ## Iteratively reveal hypotheses with the order determined by GAM
  for (i in 1:length(unrevealed_id)){
      mdl <- gam(revealed_sign ~ ns(U, df) + abs(W),family = "binomial")
      fitted.pval <- mdl$fitted.values
      fitted.pval <- fitted.pval[unrevealed_id]
      predicted.sign <- fitted.pval>0.5

      ## Reveal the W_j with smallest probability of being a positive
      ind.min <- which(fitted.pval == min(fitted.pval))
      if(length(ind.min)==1){
        ind.reveal <- ind.min
      }else{
        ind.reveal <- ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal <- unrevealed_id[ind.reveal]
      revealed_id <- c(revealed_id,ind.reveal)
      unrevealed_id <- all_id[-revealed_id]
      revealed_sign[ind.reveal] <- W_sign[ind.reveal]
      tau[ind.reveal] <- W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W, tau, offset = offset)
      if(mute == "FALSE"){print(fdphat)}
      if(fdphat <= fdr){break}
    }
    rej = which(W >= tau)

  return(rej)
}


filter_gam_variant1 <- function(W, U, U0, fdr = 0.1, offset = 1,
                                df = 5, reveal_prop = 0.5, mute = TRUE){

  ## Determine dimensionality
  p <- length(W)
  
  ## Initilization
  tau <- rep(0,p)
  rejs <- list()
  W_abs <- abs(W)
  W_sign <- as.numeric(W>0)
  revealed_sign <- rep(1,p)
  all_id <- 1:p
  tau.sel <- c()

  ## Reveal the features with U_j >= U_0
  revealed_id = which(U >= U0 | abs(W) <= quantile(abs(W), reveal_prop))

  ## Update the revealed information
  revealed_sign[revealed_id] <- W_sign[revealed_id]
  unrevealed_id <- all_id[-revealed_id]
  tau[revealed_id] <- W_abs[revealed_id] + 1

  ## Iteratively reveal hypotheses with the order determined by GAM
  for (i in 1:length(unrevealed_id)){
      mdl <- gam(revealed_sign ~ ns(U, df) + abs(W),family = "binomial")
      fitted.pval <- mdl$fitted.values
      fitted.pval <- fitted.pval[unrevealed_id]
      predicted.sign <- fitted.pval>0.5

      ## Reveal the W_j with smallest probability of being a positive
      ind.min <- which(fitted.pval == min(fitted.pval))
      if(length(ind.min)==1){
        ind.reveal <- ind.min
      }else{
        ind.reveal <- ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal <- unrevealed_id[ind.reveal]
      revealed_id <- c(revealed_id,ind.reveal)
      unrevealed_id <- all_id[-revealed_id]
      revealed_sign[ind.reveal] <- W_sign[ind.reveal]
      tau[ind.reveal] <- W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W, tau, offset = offset)
      if(mute == "FALSE"){print(fdphat)}
      if(fdphat <= fdr){break}
    }
    rej = which(W >= tau)

  return(rej)
}


calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/sum(W>=tau)
  return(fdphat)
}


## Ordering determined by a combination of loss function and side information
filter_gam_variant2 <- function(W, U, U0, lambda, fdr = 0.1, offset = 1,
                                df = 5, reveal_prop = 0.5, mute = TRUE){

  ## Determine dimensionality
  p <- length(W)
  
  ## Initilization
  tau <- rep(0,p)
  rejs <- list()
  W_abs <- abs(W)
  W_sign <- as.numeric(W>0)
  revealed_sign <- rep(1,p)
  all_id <- 1:p
  tau.sel <- c()

  ## Reveal the features with U_j >= U_0
  revealed_id = which(abs(W) <= quantile(abs(W), reveal_prop))

  ## Update the revealed information
  revealed_sign[revealed_id] <- W_sign[revealed_id]
  unrevealed_id <- all_id[-revealed_id]
  tau[revealed_id] <- W_abs[revealed_id] + 1

  ## Iteratively reveal hypotheses with the order determined by GAM
  for (i in 1:length(unrevealed_id)){
      mdl <- gam(revealed_sign ~ ns(U, df) + abs(W),family = "binomial")
      fitted.pval <- mdl$fitted.values
      fitted.pval <- (1-lambda) * fitted.pval[unrevealed_id] + 
        lambda * (U[unrevealed_id] <= U0) / 100

      ## Reveal the W_j with smallest probability of being a positive
      ind.min <- which(fitted.pval == min(fitted.pval))
      if(length(ind.min)==1){
        ind.reveal <- ind.min
      }else{
        ind.reveal <- ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal <- unrevealed_id[ind.reveal]
      revealed_id <- c(revealed_id,ind.reveal)
      unrevealed_id <- all_id[-revealed_id]
      revealed_sign[ind.reveal] <- W_sign[ind.reveal]
      tau[ind.reveal] <- W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W, tau, offset = offset)
      if(mute == "FALSE"){print(fdphat)}
      if(fdphat <= fdr){break}
    }
    rej = which(W >= tau)

  return(rej)
}


calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/sum(W>=tau)
  return(fdphat)
}

