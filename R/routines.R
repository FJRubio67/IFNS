#########################
# Some hazard functions
#########################

#--------------------------------------------------------------------------------------------------------------------------
#' Weibull Hazard
#' @param t: positive value
#' @param scale: scale parameter
#' @param shape: shape parameter
#' @return Weibull hazard function.
#' @export
hw <- function(t,scale,shape){
  pdf0 <-  dweibull(t,scale=scale,shape=shape)
  s0 <- pweibull(t,scale=scale,shape=shape, lower.tail = FALSE)
  return(pdf0/s0)
}                                                                                      

#' Weibull Cumulative hazard
#' @param t: positive value
#' @param scale: scale parameter
#' @param shape: shape parameter
#' @return Weibull cumulative hazard function.
#' @export
chw <- function(t,scale,shape){
  H0 <- -pweibull(t,scale=scale,shape=shape, lower.tail = FALSE, log.p = TRUE)  
  return(H0)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Lognormal Hazard
#' @param t: positive parameter
#' @param mu: log-location parameter
#' @param  sigma: scale parameter
#' @param  log: return log harzard (TRUE or FALSE)
#' @return lognormal hazard function
#' @export
hlnorm <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlnorm(t,mu,sigma, log = T)
  ls0 <- plnorm(t,mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}                                                                                      

#' Lognormal Cumulative hazard
#' @param t: positive parameter
#' @param mu: log-location parameter
#' @param  sigma: scale parameter
#' @return lognormal cumulative hazard function
#' @export
chlnorm <- function(t,mu,sigma){
  H0 <- -plnorm(t,mu,sigma, lower.tail = FALSE, log.p = TRUE)  
  return(H0)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Gamma Hazard
#' @param t: positive parameter
#' @param shape: shape parameter
#' @param  scale: scale parameter
#' @param  log: return log harzard (TRUE or FALSE)
#' @return Gamma hazard function
#' @export
hgamma <- function(t, shape, scale, log = FALSE){
  lpdf0 <-  dgamma(t, shape = shape, scale = scale, log = T)
  ls0 <- pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}                                                                                      

#' Gamma Cumulative hazard
#' @param t: positive parameter
#' @param shape: shape parameter
#' @param  scale: scale parameter
#' @param  log: return log harzard (TRUE or FALSE)
#' @return Gamma cumulative hazard function
#' @export
chgamma <- function(t, shape, scale){
  H0 <- -pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)  
  return(H0)
}

#--------------------------------------------------------------------------------------------------------------------------
#' PGW Survival Function
#' @param t: positive parameter
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @param  log.p: return log survival (TRUE or FALSE)
#' @return PGW survival function
#' @export
spgw <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- 1 - ( 1 + (t/sigma)^nu )^(1/gamma)
  if(log.p) return(val) else return(exp(val))
}

#' PGW Hazard Function
#' @param t: positive parameter
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @param  log: return log hazard (TRUE or FALSE)
#' @return PGW hazard function
#' @export
hpgw <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(nu) - log(gamma) - nu*log(sigma) + (nu-1)*log(t) + 
    (1/gamma - 1)*log( 1 + (t/sigma)^nu )
  if(log) return(val) else return(exp(val))
}

#' PGW Cumulative Hazard Function
#' @param t: positive parameter
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @return PGW cumulative hazard function
#' @export
chpgw <- function(t, sigma, nu, gamma){
  val <- -1 + ( 1 + (t/sigma)^nu )^(1/gamma)
  return(val) 
}

#' Quantile Function
#' @param p: probability value
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @return vector of quantiles
#' @export
qpgw <- function(p, sigma, nu, gamma){
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(out)
}

#' Random Number Generation Function PGW distribution
#' @param n: sample size
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @return vector of simulated values
#' @export
rpgw <- function(n, sigma, nu, gamma){
  p <- runif(n)
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(as.vector(out))
}

###########################################################################################################
#' Function to simulate times to event from a model with PGW baseline hazard
#' and GH structure
###########################################################################################################
#' @param seed: simulation seed
#' @param n: number of simulations
#' @param des_t: design matrix for time-level effects
#' @param des: design matrix for hazard-level effects
#' @param parPGW: parameters of the PGW baseline hazard
#' @param alpha: regression coefficients for time-level effects
#' @param beta: regression coefficients for hazard-level effects
#' @param frailty: "One" indicates frailties are one (no frailty), "Gamma" indicates gamma-distributed frailties
#' @param scale.frail: scale parameter for the frailty distribution
#' @return returns the vector "lifestimes" containing the simulated survival times. 
#' @export 
simPGWGH <- function(seed, n, des_t, des, parPGW, alpha, beta, frailty = "One", scale.frail = NULL){
  a0 <- parPGW[1]; b0 <- parPGW[2]; c0 <- parPGW[3];
  des_t <- as.matrix(des_t); des <- as.matrix(des)
  X.alpha <- des_t%*%alpha
  S.beta <- S%*%beta
  exp.xalpha  <- as.vector(exp(X.alpha))
  exp.dif <- as.vector(exp(X.alpha-S.beta))
  # Simulation
  set.seed(seed)
  u  <- runif(n)
  if(frailty=="One"){
    val <- 1-exp( log( 1-u )*exp.dif )
    lifestimes <- qpgw(val , a0, b0, c0)/exp.xalpha
  }
  if(frailty=="Gamma"){
    lambda <- rgamma(n, scale = scale.frail, shape = 1/scale.frail)
    lifestimes <-  val <- 1-exp( log( 1-u )*exp.dif/lambda )
    lifestimes <- qpgw(val , a0, b0, c0)/exp.xalpha
  }
  return(as.vector(lifestimes))
}

###########################################################################################################
#' Function to simulate times to event from a model with LN baseline hazard
#' and GH structure
###########################################################################################################
#' @param seed: simulation seed
#' @param n: number of simulations
#' @param des_t: design matrix for time-level effects
#' @param des: design matrix for hazard-level effects
#' @param parLN: parameters of the lognormal baseline hazard
#' @param alpha: regression coefficients for time-level effects
#' @param beta: regression coefficients for hazard-level effects
#' @param frailty: "One" indicates frailties are one (no frailty), "Gamma" indicates gamma-distributed frailties
#' @param scale.frail: scale parameter for the frailty distribution
#' @return returns the vector "lifestimes" containing the simulated survival times. 
#' @export 
simLNGH <- function(seed, n, X, S, parLN, alpha, beta, frailty = "One", scale.frail = NULL){
  a0 <- parLN[1]; b0 <- parLN[2]; 
  X <- as.matrix(X); S <- as.matrix(S)
  X.alpha <- X%*%alpha
  S.beta <- S%*%beta
  exp.xalpha  <- as.vector(exp(X.alpha))
  exp.dif <- as.vector(exp(X.alpha-S.beta))
  # Simulation
  set.seed(seed)
  u  <- runif(n)
  if(frailty=="One"){
    val <- 1-exp( log( 1-u )*exp.dif )
    lifestimes <- qlnorm(val , a0, b0)/exp.xalpha
  }
  if(frailty=="Gamma"){
    lambda <- rgamma(n, scale = scale.frail, shape = 1/scale.frail)
    lifestimes <-  val <- 1-exp( log( 1-u )*exp.dif/lambda )
    lifestimes <- qlnorm(val , a0, b0)/exp.xalpha
  }
  return(as.vector(lifestimes))
}

#----------------------------------------------------------------------
#' Excess Hazard Models: Weibull baseline
#' WMLE function. It returns the output from optim or nlminb
#----------------------------------------------------------------------

#' @param init    : initial point for optimisation step   (log(scale), log(shape), beta), 
#' @param hstr    : hazard structure, Weibull without covariates ("W")
#' Weibull with AFT model (WAFT)
#' Weibull with PH model (WPH)
#' Weibull with AH model (WAH)
#' @param method  : "nlminb" or optimisation method to be used in optim 
#' @param maxit   : maximum number of iterations in optim or nlminb
#' @param times   : times to event
#' @param status    : status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param rates   : expected population mortality hazard rates
#' @param des     : design matrix  
#' @return WMLE returns an object containing the output of optim or nlminb
#' @export
WMLE <- function(init, hstr = "W", times, status, rates, des = NULL, method = "Nelder-Mead", maxit = 100){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des))  des <- as.matrix(des); des.obs <- des[status,]
  
  # Weibull model
  if(hstr == "W"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);
      haz0 <- hp.as.obs + hw(times.obs,ae0,be0)
      val <- - sum(log(haz0)) + sum(chw(times,ae0,be0))
      return(val)
    }
  } 
  # PH model  
  if(hstr == "WPH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]); beta0 <- par[3:(3+p-1)];
      exp.x.beta <- as.vector(exp(des%*%beta0))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hw(times.obs,ae0,be0)*exp.x.beta.obs
      val <- -sum(log(haz0)) + sum(chw(times,ae0,be0)*exp.x.beta)
      return(val)
    }
  } 
  # AFT Model  
  if(hstr == "WAFT"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]); beta0 <- par[3:(3+p-1)];
      exp.x.beta <- as.vector(exp(des%*%beta0))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hw(times.obs*exp.x.beta.obs,ae0,be0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chw(times*exp.x.beta,ae0,be0))
      return(val)
    }
  }
  if(hstr == "WAH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]); beta0 <- par[3:(3+p-1)];
      exp.x.beta <- as.vector(exp(des%*%beta0))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hw(times.obs*exp.x.beta.obs,ae0,be0)
      val <- - sum(log(haz0)) + sum(chw(times*exp.x.beta,ae0,be0)/exp.x.beta)
      return(sum(val))
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit),method=method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  return(OPT)
}

#----------------------------------------------------------------------
#' Excess Hazard Models: PGW baseline with GH structure
#' PGWMLE function. It returns the output from optim or nlminb
#----------------------------------------------------------------------

#' @param init    : initial point for optimisation step
#'            (log(scale), log(shape1), log(shape2), alpha, beta), 
#'                       where alpha is only required for PGWGH   
#' @param hstr    : hazard structure, PGW without covariates ("PGW"),
#'            PGW with AFT model (PGWAFT)
#'            PGW with PH model (PGWPH)
#'            PGW with AH model (PGWAH)
#'            PGW with GH model (PGWGH)
#' @param method  : "nlminb" or optimisation method to be used in optim  
#' @param maxit   : maximum number of iterations in optim
#' @param times   : times to event
#' @param status    : status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param rates   : expected population mortality hazard rates
#' @param des     : design matrix for proportional hazard effects 
#' @param des_t   : design matrix for time-dependent effects (it is recommended not to use splines here)
#' @return list containing the output from the optimsation (OPT) and the log-likelihood (loglik)
#' @export

PGWMLE <- function(init, times, status, rates, hstr = "PGW", des = NULL, des_t = NULL, method = "Nelder-Mead", maxit = 100){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  
  # PGW model
  if(hstr == "PGW"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]);
      haz0 <- hp.as.obs + hpgw(times.obs,ae0,be0,ce0)
      val <- - sum(log(haz0)) + sum(chpgw(times,ae0,be0,ce0))
      return(val)
    }
  } 
  # PH model  
  if(hstr == "PGWPH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); beta <- par[4:(3+p)];
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs,ae0,be0,ce0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chpgw(times,ae0,be0,ce0)*exp.x.beta)
      return(val)
    }
  } 
  # AFT Model  
  if(hstr == "PGWAFT"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); beta <- par[4:(3+p)];
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.beta.obs,ae0,be0,ce0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chpgw(times*exp.x.beta,ae0,be0,ce0))
      return(val)
    }
  }
  # AH Model  
  if(hstr == "PGWAH"){
    p <- ncol(des_t)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p)];
      exp.x.alpha <- as.vector(exp(des_t%*%alpha))
      exp.x.alpha.obs <- exp.x.alpha[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)
      val <- - sum(log(haz0)) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)/exp.x.alpha)
      return(val)
    }
  }
  # GH Model  
  if(hstr == "PGWGH"){
    p0 <- dim(des_t)[2]
    p1 <- dim(des)[2]
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      exp.x.alpha <- as.vector(exp(des_t%*%alpha))
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.dif <- as.vector(exp( des%*%beta - des_t%*%alpha ))
      exp.x.alpha.obs <- exp.x.alpha[status]
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit),method=method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  out <- list(OPT = OPT, loglik = log.lik)
  return(out)
}

#----------------------------------------------------------------------
#' Excess Hazard Models: LogNormal baseline with GH structure
#' LNMLE function. It returns the output from optim or nlminb
#----------------------------------------------------------------------

#' @param init    : initial point for optimisation step
#'           (log-location, log(scale), alpha, beta), 
#'           where alpha is only required for LNGH   
#' @param hstr    : hazard structure, LN without covariates ("LN"),
#'           LN with AFT model (LNAFT)
#'           LN with PH model (LNPH)
#'           LN with AH model (LNAH)
#'           LN with GH model (LNGH)
#' @param method  : "nlminb" or optimisation method to be used in optim 
#' @param maxit   : maximum number of iterations in optim
#' @param times   : times to event
#' @param status    : status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param rates   : expected population mortality hazard rates
#' @param des     : design matrix for proportional hazard effects 
#' @param des_t   : design matrix for time-dependent effects (it is recommended not to use splines here)
#' @return list containing the output from the optimsation (OPT) and the log-likelihood (loglik)
#' @export

LNMLE <- function(init, times, status, rates, hstr = "LN", des = NULL, des_t = NULL, method = "Nelder-Mead", maxit = 100){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  
  # LN model
  if(hstr == "LN"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);
      haz0 <- hp.as.obs + hlnorm(times.obs,ae0,be0)
      val <- - sum(log(haz0)) + sum(chlnorm(times,ae0,be0))
      return(val)
    }
  } 
  # PH model  
  if(hstr == "LNPH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);  beta <- par[3:(2+p)];
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs,ae0,be0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chlnorm(times,ae0,be0)*exp.x.beta)
      return(val)
    }
  }
  # AFT Model  
  if(hstr == "LNAFT"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);  beta <- par[3:(2+p)];
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.beta.obs,ae0,be0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chlnorm(times*exp.x.beta,ae0,be0))
      return(val)
    }
  }
  # AH Model  
  if(hstr == "LNAH"){
    p <- ncol(des_t)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);  alpha <- par[3:(2+p)];
      exp.x.beta <- as.vector(exp(des_t%*%alpha))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.beta.obs,ae0,be0)
      val <- - sum(log(haz0)) + sum(chlnorm(times*exp.x.beta,ae0,be0)/exp.x.beta)
      return(val)
    }
  }
  # GH Model  
  if(hstr == "LNGH"){
    p0 <- dim(des_t)[2]
    p1 <- dim(des)[2]
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);  alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      exp.x.alpha <- as.vector(exp(des_t%*%alpha))
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.dif <- as.vector(exp( des%*%beta - des_t%*%alpha ))
      exp.x.alpha.obs <- exp.x.alpha[status]
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.alpha.obs,ae0,be0)*exp.x.beta.obs
      val <- - sum(log(haz0)) + sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit),method=method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  out <- list(OPT = OPT, loglik = log.lik)
  return(out)
}

#----------------------------------------------------------------------
#' Individual Fraitly Excess Hazard Models: PGW Baseline hazard
#' IFPGWMLE function. It returns the output from optim or nlminb
#----------------------------------------------------------------------

#' @param init   : initial point for optimisation step
#'           (log(scale), log(shape1), log(shape2), alpha, beta, log(b)), 
#'           where alpha is only required for PGWGH   
#' @param hstr   : hazard structure, PGW without covariates ("IFPGW"),
#'          PGW with AFT model (IFPGWAFT)
#'          PGW with PH model (IFPGWPH)
#'          PGW with AH model (IFPGWAH)
#'          PGW with GH model (IFPGWGH)
#' @param times   : times to event
#' @param status    : status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param rates   : expected population mortality hazard rates
#' @param method : optimisation method to be used in optim 
#' @param maxit  : maximum number of iterations in optim
#' @return list containing the output from the optimsation (OPT) and the log-likelihood (loglik)
#' @export

IFPGWMLE <- function(init, hstr = "IFPGW", times, status, rates, des = NULL, des_t = NULL, method = "Nelder-Mead", maxit = 100){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  
  # IFPGW model
  if(hstr == "IFPGW"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]);  b0 <- exp(par[4])
      haz0 <- hp.as.obs + hpgw(times.obs,ae0,be0,ce0)/(  1 + b0*chpgw(times.obs,ae0,be0,ce0))
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times,ae0,be0,ce0) ))
      return(val)
    }
  } 
  # IFPH model  
  if(hstr == "IFPGWPH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); beta <- par[4:(3+p)];  b0 <- exp(par[p+4])
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs,ae0,be0,ce0)*exp.x.beta.obs/(  1 + b0*chpgw(times.obs,ae0,be0,ce0)*exp.x.beta.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times,ae0,be0,ce0)*exp.x.beta))
      return(val)
    }
  } 
  # IFAFT Model  
  if(hstr == "IFPGWAFT"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); beta <- par[4:(3+p)];  b0 <- exp(par[p+4])
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.beta.obs,ae0,be0,ce0)*exp.x.beta.obs/(  1 + b0*chpgw(times.obs*exp.x.beta.obs,ae0,be0,ce0))
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times*exp.x.beta,ae0,be0,ce0)))
      return(val)
    }
  }
  # IFAH Model
  if(hstr == "IFPGWAH"){
    p <- ncol(des_t)
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p)];  b0 <- exp(par[p+4])
      exp.x.beta <- as.vector(exp(des_t%*%alpha))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.beta.obs,ae0,be0,ce0)/(  1 + b0*chpgw(times.obs*exp.x.beta.obs,ae0,be0,ce0)/exp.x.beta.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times*exp.x.beta,ae0,be0,ce0)/exp.x.beta))
      return(val)
    }
  }
  # IFGH Model
  if(hstr == "IFPGWGH"){
    p0 <- dim(des_t)[2]
    p1 <- dim(des)[2]
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)];  b0 <- exp(par[p0+p1+4])
      x.alpha <- des_t%*%alpha
      x.beta <- des%*%beta
      exp.x.alpha <- as.vector(exp(x.alpha))
      exp.x.beta <- as.vector(exp(x.beta ))
      exp.x.beta.dif <- as.vector(exp(x.beta-x.alpha))
      exp.x.alpha.obs <- exp.x.alpha[status]
      exp.x.beta.obs <- exp.x.beta[status]
      exp.x.beta.dif.obs <- exp.x.beta.dif[status]
      haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs/(  1 + b0*chpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.dif.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif))
      return(val)
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit),method=method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  out <- list(OPT = OPT, loglik = log.lik)
  return(out)
}

#----------------------------------------------------------------------
#' Individual Fraitly Excess Hazard Models: LN Baseline hazard
#' IFLNMLE function. It returns the output from optim or nlminb
#----------------------------------------------------------------------

#' @param init   : initial point for optimisation step
#'           (log-location, log(scale), alpha, beta, log(b)), 
#'           where alpha is only required for LNGH   
#' @param hstr   : hazard structure, LN without covariates ("IFLN"),
#'          LN with AFT model (IFLNAFT)
#'          LN with PH model (IFLNPH)
#'          LN with AH model (IFLNAH)
#'          LN with GH model (IFLNGH)
#' @param times   : times to event
#' @param status    : status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param rates   : expected population mortality hazard rates
#' @param method : optimisation method to be used in optim 
#' @param maxit  : maximum number of iterations in optim
#' @return list containing the output from the optimsation (OPT) and the log-likelihood (loglik)
#' @export

IFLNMLE <- function(init, hstr = "IFLN", times, status, rates, des = NULL, des_t = NULL, method = "Nelder-Mead", maxit = 100){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  
  # IFLN model
  if(hstr == "IFLN"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);  b0 <- exp(par[3])
      haz0 <- hp.as.obs + hlnorm(times.obs,ae0,be0)/(  1 + b0*chlnorm(times.obs,ae0,be0))
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chlnorm(times,ae0,be0) ))
      return(val)
    }
  } 
  # IFPH model  
  if(hstr == "IFLNPH"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   beta <- par[3:(2+p)];  b0 <- exp(par[p+3])
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs,ae0,be0)*exp.x.beta.obs/(  1 + b0*chlnorm(times.obs,ae0,be0)*exp.x.beta.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chlnorm(times,ae0,be0)*exp.x.beta))
      return(val)
    }
  } 
  # IFAFT Model  
  if(hstr == "IFLNAFT"){
    p <- ncol(des)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   beta <- par[3:(2+p)];  b0 <- exp(par[p+3])
      exp.x.beta <- as.vector(exp(des%*%beta))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.beta.obs,ae0,be0)*exp.x.beta.obs/(  1 + b0*chlnorm(times.obs*exp.x.beta.obs,ae0,be0))
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chlnorm(times*exp.x.beta,ae0,be0)))
      return(val)
    }
  }
  # IFAH Model
  if(hstr == "IFLNAH"){
    p <- ncol(des_t)
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p)];  b0 <- exp(par[p+3])
      exp.x.beta <- as.vector(exp(des_t%*%alpha))
      exp.x.beta.obs <- exp.x.beta[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.beta.obs,ae0,be0)/(  1 + b0*chlnorm(times.obs*exp.x.beta.obs,ae0,be0)/exp.x.beta.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chlnorm(times*exp.x.beta,ae0,be0)/exp.x.beta))
      return(val)
    }
  }
  # IFGH Model
  if(hstr == "IFLNGH"){
    p0 <- dim(des_t)[2]
    p1 <- dim(des)[2]
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)];  b0 <- exp(par[p0+p1+3])
      x.alpha <- des_t%*%alpha
      x.beta <- des%*%beta
      exp.x.alpha <- as.vector(exp(x.alpha))
      exp.x.beta <- as.vector(exp(x.beta ))
      exp.x.beta.dif <- as.vector(exp(x.beta-x.alpha))
      exp.x.alpha.obs <- exp.x.alpha[status]
      exp.x.beta.obs <- exp.x.beta[status]
      exp.x.beta.dif.obs <- exp.x.beta.dif[status]
      haz0 <- hp.as.obs + hlnorm(times.obs*exp.x.alpha.obs,ae0,be0)*exp.x.beta.obs/(  1 + b0*chlnorm(times.obs*exp.x.alpha.obs,ae0,be0)*exp.x.beta.dif.obs)
      val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif))
      return(val)
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit),method=method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  out <- list(OPT = OPT, loglik = log.lik)
  return(out)
}

###########################################################################################
#' Function to calculate the normal confidence intervals on the log scale for positive par
###########################################################################################
#' @param FUN   : minus log-likelihood function to be used to calculate the confidence intervals
#' @param MLE   : maximum likelihood estimator of the parameters of interest
#' @param level : confidence level
#' @return matrix with colums ("Lower","Upper","Transf MLE", "Std. Error")
#' @export

Conf.Int <- function(FUN,MLE,level=0.95){
  sd.int <- abs(qnorm(0.5*(1-level)))
  HESS <- hessian(FUN,x=MLE)
  Fisher.Info <- solve(HESS)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- MLE + sd.int*Sigma
  L<- MLE - sd.int*Sigma
  C.I <- cbind(L,U,MLE, Sigma)
  rownames(C.I)<- paste0("par", seq_along(1:length(MLE)))
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

###########################################################################################
#' Function to calculate std errors and inclusion of the true parameter values
###########################################################################################
#' @param trans.mle : mle with log scale for positive parameters
#' @param level: confidence level
#' @param times: survival times
#' @param status: vital status
#' @param rates: population hazard values
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return confidence intervals for PGWGH model
#' @export
Conf.IntPGWGH <- function(trans.mle, level = 0.95, times, status, rates, des = NULL, des_t = NULL){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des_t))  des <- as.matrix(des); des.obs <- des[status,];
  if(!is.null(des_t))  des_t <- as.matrix(des_t); des_t.obs <- des_t[status,];
  # GH Model  
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]
  log.lik <- function(par){
    ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
    exp.x.alpha <- as.vector(exp(des_t%*%alpha))
    exp.x.beta <- as.vector(exp(des%*%beta))
    exp.x.beta.dif <- as.vector(exp( des%*%beta - des_t%*%alpha ))
    exp.x.alpha.obs <- exp.x.alpha[status]
    exp.x.beta.obs <- exp.x.beta[status]
    haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs
    val <- - sum(log(haz0)) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
    return(sum(val))
  }
  CI <- Conf.Int(log.lik,trans.mle,level=0.95)
  CI[1:3,1:2] <- exp(CI[1:3,1:2])
  # Standard Error and coverage 
  sd.error <- CI[,4]
  sd.error[1:3] <- sqrt(exp(2*CI[1:3,3])*sd.error[1:3]^2)
  counter <- vector()
  for(i in 1:(p0+p1+3)) counter[i] <- ifelse( true.par[i] < CI[i,2] & true.par[i] > CI[i,1],  1,0)
  out <- as.matrix(cbind(sd.error,counter))
  colnames(out) <- c("Std. Error", "Coverage")
  return(out)
}

####################################################################################################################
#' Function to calculate std errors and inclusion of the true parameter values: frailty model
####################################################################################################################
#' @param trans.mle : mle with log scale for positive parameters
#' @param level: confidence level
#' @param times: survival times
#' @param status: vital status
#' @param rates: population hazard values
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return confidence intervals for IFPGWGH model
#' @export
IFConf.IntPGWGH <- function(trans.mle, level = 0.95, times, status, rates, des = NULL, des_t = NULL){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des))  des <- as.matrix(des); des.obs <- des[status,];
  if(!is.null(des_t))  des_t <- as.matrix(des_t); des_t.obs <- des_t[status,];
  # GH Model  
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]
  log.lik <- function(par){
    ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)];  b0 <- exp(par[p0+p1+4])
    x.alpha <- des_t%*%alpha
    x.beta <- des%*%beta
    exp.x.alpha <- as.vector(exp(x.alpha))
    exp.x.beta <- as.vector(exp(x.beta ))
    exp.x.beta.dif <- as.vector(exp(x.beta-x.alpha))
    exp.x.alpha.obs <- exp.x.alpha[status]
    exp.x.beta.obs <- exp.x.beta[status]
    exp.x.beta.dif.obs <- exp.x.beta.dif[status]
    haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs/(  1 + b0*chpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.dif.obs)
    val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif))
    return(val)
  }
  CI <- Conf.Int(log.lik,trans.mle,level=0.95)
  CI[1:3,1:2] <- exp(CI[1:3,1:2])
  CI[p0+p1+4,1:2] <- exp(CI[p0+p1+4,1:2])
  # Standard Error and coverage 
  sd.error <- CI[,4]
  sd.error[1:3] <- sqrt(exp(2*CI[1:3,3])*sd.error[1:3]^2)
  sd.error[p0+p1+4] <- sqrt(exp(2*CI[p0+p1+4,3])*sd.error[p0+p1+4]^2)
  counter <- vector()
  for(i in 1:(p0+p1+4)) counter[i] <- ifelse( true.parIF[i] < CI[i,2] & true.parIF[i] > CI[i,1],  1,0)
  out <- as.matrix(cbind(sd.error,counter))
  colnames(out) <- c("Std. Error", "Coverage")
  return(out)
}

#' Function to evaluate the net survival based on specific values of the parameters and covariates
#' @param t: positive value
#' @param sigma: scale parameter
#' @param  nu: shape parameter
#' @param gamma: shape parameter
#' @param alpha: regression parameters for time-level effects
#' @param beta: regression coefficients for hazard-level effects
#' @param scale.frail: scale parameter of the frailty distribution
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return Net survival for IFPGW model: IFPGW, IFPGWPH, IFPGWAH, IFPGWAFT, IFPGWGH
#' @export

IFNS <- function(t, sigma, nu, gamma, alpha = NULL, beta = NULL, scale.frail, des=NULL, des_t=NULL, hstr = "IFPGW"){
  # IFPGW model
  if(hstr == "IFPGW"){
    s.val <-   1/( 1 + scale.frail*chpgw(t, sigma, nu, gamma) )^(1/scale.frail)
  }
  # IFPGW-PH model
  if(hstr == "IFPGWPH"){
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t, sigma, nu, gamma)*exp(x.beta) )^(1/scale.frail)
  }
  # IFPGW-AH model
  if(hstr == "IFPGWAH"){
    x.alpha <- des_t%*%alpha
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.alpha), sigma, nu, gamma)*exp(- x.alpha) )^(1/scale.frail)
  }
  # IFPGW-AFT model
  if(hstr == "IFPGWAFT"){
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.beta), sigma, nu, gamma) )^(1/scale.frail)
  }
  # IFPGW-GH model
  if(hstr == "IFPGWGH"){
    x.alpha <- des_t%*%alpha
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.alpha), sigma, nu, gamma)*exp(x.beta - x.alpha) )^(1/scale.frail)
  }
  return(mean(s.val))
}


# Function to evaluate the net survival based on specific values of the parameters and covariates for each individual
IFNS.ind <- function(t, sigma, nu, gamma, alpha = NULL, beta = NULL, scale.frail, des=NULL, des_t=NULL, hstr = "IFPGW"){
  if(!is.null(des))  des <- as.matrix(des); des.obs <- des[status,];
  if(!is.null(des_t))  des_t <- as.matrix(des_t); des_t.obs <- des_t[status,];
  # IFPGW model
  if(hstr == "IFPGW"){
    s.val <-   1/( 1 + scale.frail*chpgw(t, sigma, nu, gamma) )^(1/scale.frail)
  }
  # IFPGW-PH model
  if(hstr == "IFPGWPH"){
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t, sigma, nu, gamma)*exp(x.beta) )^(1/scale.frail)
  }
  # IFPGW-AH model
  if(hstr == "IFPGWAH"){
    x.alpha <- des_t%*%alpha
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.alpha), sigma, nu, gamma)*exp(- x.alpha) )^(1/scale.frail)
  }
  # IFPGW-AFT model
  if(hstr == "IFPGWAFT"){
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.beta), sigma, nu, gamma) )^(1/scale.frail)
  }
  # IFPGW-GH model
  if(hstr == "IFPGWGH"){
    x.alpha <- des_t%*%alpha
    x.beta <- des%*%beta
    s.val <-   1/( 1 + scale.frail*chpgw(t*exp(x.alpha), sigma, nu, gamma)*exp(x.beta - x.alpha) )^(1/scale.frail)
  }
  return(s.val)
}




######################################################
#' Net survival function for PGWGH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return net survival for PGWGH model
#' @export
NSfit_PGWGH <- function(t,par,des_t,des){
  des_t <- as.matrix(des_t); des <- as.matrix(des)
  # Classical model
  exp.mle0 <- as.vector(exp( des_t%*%par[4:(3+ncol(des_t))] ))
  exp.mledif0 <- as.vector(exp( des%*%par[(4+ncol(des_t)):(3+ncol(des_t)+ncol(des))]  -des_t%*%par[4:(3+ncol(des_t))] ))
  out <- mean(exp(-chpgw(t*exp.mle0,par[1],par[2],par[3])*exp.mledif0))  
  return(out)
}


######################################################
#' Net survival function for IFPGWGH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return net survival for IFPGWGH model
#' @export
NSfit_IFPGWGH <- function(t,par,des_t,des){
  des_t <- as.matrix(des_t); des <- as.matrix(des)
  # Frailty model
  exp.mle <- as.vector(exp( des_t%*%par[4:(3+ncol(des_t))] ))
  exp.mledif <- as.vector(exp( des%*%par[(4+ncol(des_t)):(3+ncol(des_t)+ncol(des))] - des_t%*%par[4:(3+ncol(des_t))] ))
  bfr <- tail(par, n=1)
  out <-  mean(1/(1+bfr*chpgw(t*exp.mle,par[1],par[2],par[3])*exp.mledif)^(1/bfr) )
  return(out)
}


######################################################
#' Net survival function for PGWPH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @return net survival for PGWPH model
#' @export
NSfit_PGWPH <- function(t,par,des){
des <- as.matrix(des)
  # Classical model
  exp.mle <- as.vector(exp( des%*%par[(4):(3+ncol(des))] ))
  out <- mean(exp(-chpgw(t,par[1],par[2],par[3])*exp.mle))  
  return(out)
}


######################################################
#' Net survival function for IFPGWPH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @return net survival for IFPGWPH model
#' @export
NSfit_IFPGWPH <- function(t,par,des){
  des <- as.matrix(des)
  # Frailty model
  exp.mle <- as.vector(exp( des%*%par[(4):(3+ncol(des))]))
  bfr <- tail(par, n=1)
  out <-  mean(1/(1+bfr*chpgw(t,par[1],par[2],par[3])*exp.mle)^(1/bfr) )
  return(out)
}


######################################################
#' Net survival function for LNGH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return net survival for LNGH model
#' @export
NSfit_LNGH <- function(t,par,des_t,des){
  des_t <- as.matrix(des_t); des <- as.matrix(des)
  # Classical model
  exp.mle0 <- as.vector(exp( des_t%*%par[3:(2+ncol(des_t))] ))
  exp.mledif0 <- as.vector(exp( des%*%par[(3+ncol(des_t)):(2+ncol(des_t)+ncol(des))]  -des_t%*%par[3:(2+ncol(des_t))] ))
  out <- mean(exp(-chlnorm(t*exp.mle0,par[1],par[2])*exp.mledif0))  
  return(out)
}


######################################################
#' Net survival function for IFLNGH model
######################################################
#' @param t: positive value
#' @param par: model parameters
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @return net survival for IFLNGH model
#' @export
NSfit_IFLNGH <- function(t,par,des_t,des){
  des_t <- as.matrix(des_t); des <- as.matrix(des)
  # Frailty model
  exp.mle <- as.vector(exp( des_t%*%par[3:(2+ncol(des_t))] ))
  exp.mledif <- as.vector(exp( des%*%par[(3+ncol(des_t)):(2+ncol(des_t)+ncol(des))] - des_t%*%par[3:(2+ncol(des_t))] ))
  bfr <- tail(par, n=1)
  out <-  mean(1/(1+bfr*chlnorm(t*exp.mle,par[1],par[2])*exp.mledif)^(1/bfr) )
  return(out)
}


#########################################################################################################################
#' Confidence intervals for Net survival at time t0: PGWGH model
#########################################################################################################################

#' @param  t0 : time at what the confidence interval will be calculated
#' @param NMC : number of Monte Carlo iterations
#' @param MLE: Maximum likelihood estimate
#' @param level : confidence level
#' @param times: survival times
#' @param status: vital status
#' @param rates: population hazard values
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @param subgroup: indices of subgroup of interest
#' @return matrix containing ("lower","estimate","upper")
#' @export

CI.NSfit_PGWGH <- function(t0, NMC, MLE, level, times, rates, status, des_t, des, subgroup){
  # Required variables
  subgroup <- as.vector(subgroup)
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]
  
  # log-likelihood function
  log.lik <- function(par){
    ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
    exp.x.alpha <- as.vector(exp(des_t%*%alpha))
    exp.x.beta <- as.vector(exp(des%*%beta))
    exp.x.beta.dif <- as.vector(exp( des%*%beta - des_t%*%alpha ))
    exp.x.alpha.obs <- exp.x.alpha[status]
    exp.x.beta.obs <- exp.x.beta[status]
    haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs
    val <- - sum(log(haz0)) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
    return(sum(val))
  }
  # Covariance matrix
  r.MLE <- MLE
  r.MLE[1:3] <- log(MLE[1:3])
  HESS <- hessian(log.lik,x=r.MLE)
  Sigma <- solve(HESS)
  
  # Monte Carlo iterations
  iter <- vector()
  for(i in 1:NMC){
    val <- rmvnorm(1,mean = r.MLE, sigma = Sigma)
    val[1:3] <- exp(val[1:3])
    iter[i] <- NSfit_PGWGH(t0,val,des_t[subgroup,],des[subgroup,])
  }
  
  L <- quantile(iter,(1-level)*0.5)
  U <- quantile(iter,(1+level)*0.5)
  
  M <- NSfit_PGWGH(t0,MLE,des_t[subgroup,],des[subgroup,])
  
  out <- c(L,M,U)
  names(out) <- c("lower","estimate","upper")
  return(out)
}

#########################################################################################################################
#' Confidence intervals for Net survival at time t0: IFPGWGH model
#########################################################################################################################

#' @param  t0 : time at what the confidence interval will be calculated
#' @param NMC : number of Monte Carlo iterations
#' @param MLE: Maximum likelihood estimate
#' @param level : confidence level
#' @param times: survival times
#' @param status: vital status
#' @param rates: population hazard values
#' @param des: design matrix for hazard-level effects
#' @param des_t: design matrix for time-level effects
#' @param subgroup: indices of subgroup of interest
#' @return matrix containing ("lower","estimate","upper")
#' @export

CI.NSfit_IFPGWGH <- function(t0, NMC, MLE, level, times, rates, status, des_t, des, subgroup){
  # Required variables
  subgroup <- as.vector(subgroup)
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  hp.as <- as.vector(rates)
  hp.as.obs <- rates[status]
  times.obs <- times[status]
  if(!is.null(des)) des <- as.matrix(des); des.obs <- des[status,]
  if(!is.null(des_t)) des_t <- as.matrix(des_t); des_t.obs <- des_t[status,]
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]
  
  # log-likelihood function
  log.lik <- function(par){
    ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)];  b0 <- exp(par[p0+p1+4])
    x.alpha <- des_t%*%alpha
    x.beta <- des%*%beta
    exp.x.alpha <- as.vector(exp(x.alpha))
    exp.x.beta <- as.vector(exp(x.beta ))
    exp.x.beta.dif <- as.vector(exp(x.beta-x.alpha))
    exp.x.alpha.obs <- exp.x.alpha[status]
    exp.x.beta.obs <- exp.x.beta[status]
    exp.x.beta.dif.obs <- exp.x.beta.dif[status]
    haz0 <- hp.as.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.obs/(  1 + b0*chpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0)*exp.x.beta.dif.obs)
    val <- - sum(log(haz0)) + sum((1/b0)*log( 1 + b0*chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif))
    return(val)
  }
  # Covariance matrix
  r.MLE <- c(log(MLE[1:3]),MLE[-c(1:3,length(MLE))],log(MLE[length(MLE)]))
  
  HESS <- hessian(log.lik,x=r.MLE)
  Sigma <- solve(HESS)
  
  # Monte Carlo iterations
  iter <- vector()
  for(i in 1:NMC){
    val <- rmvnorm(1,mean = r.MLE, sigma = Sigma)
    val <- c(exp(val[1:3]),val[-c(1:3,length(val))],exp(val[length(val)]))
    iter[i] <- NSfit_IFPGWGH(t0,val,des_t[subgroup,],des[subgroup,])
  }
  
  L <- quantile(iter,(1-level)*0.5)
  U <- quantile(iter,(1+level)*0.5)
  
  M <- NSfit_IFPGWGH(t0,MLE,des_t[subgroup,],des[subgroup,])
  
  out <- c(L,M,U)
  names(out) <- c("lower","estimate","upper")
  return(out)
}


