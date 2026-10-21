#' Tumor growth model
#'
#' @param data A data frame containing the observations and variables necessary
#' for a continuous tumor growth model estimation.
#' @param variables A list containing the names of variables in the inputed data frame \code{data}.
#' \describe{
#'  \item{id}{A numeric id for each subject.}
#'  \item{det}{the cancer detection mode: 1 through symptoms, 2 through screening}
#'  \item{tdet}{the age at cancer detection}
#'  \item{ddet}{the tumor diameter (mm) at cancer detection}
#'  \item{scr}{a list of screening results, 0 for negative, 1 for positive.}
#'  \item{tscr}{a list for each patient of their age during their previous negative screening.}
#' }
#' @param d0 The smallest diameter of a tumor detectable by screening defining the starting point of the tumor growth function. Default is 0.5
#' @param distribution.t0 A str term to define the distribution of age at onset.
#' @param initpar An optional vector of p numeric terms for initializing model parameters.
#' @param fixpar An optional boolean vector with TRUE to fix model parameters
#' @param integrate.r A short term used to define the integration method.
#' “kronrod31” for a Gauss-Kronrod quadrature with 31 nodes.
#' “adaptive” for the adaptive Gauss-Kronrod quadrature from the stats package.
#' @param n_cores The number cores for parallel computation. Default to 1 (sequential mode),
#' parallelization is useful for saving time on the often long estimation of continuous tumor growth models.
#'
#' @importFrom itertools isplitVector
#' @importFrom stats dgamma dexp dweibull dlnorm dnorm pexp pgamma plnorm pweibull model.matrix na.omit integrate
#' @importFrom utils tail
#' @import doParallel
#' @import foreach
#' @import marqLevAlg
#'
#' @return A list with components:
#' \item{par}{The vector of estimated parameters.}
#' \item{parnames}{The vector containing the parameter names.}
#' \item{convergence}{Convergence status. TRUE/FALSE.}
#' \item{lld}{The value of the Log-likelihood function.}
#' \item{vcov}{The parameter variance-covariance matrix.}
#' \item{iterations}{The number of iterations of the Marquardt-Levenberg optimization}
#'
#' @export
TumorGrowthModel <- function(data,
                             variables=list(id = "id",
                                            det = "det",
                                            tdet = "tdet",
                                            ddet = "ddet",
                                            scr = "scr",
                                            tscr = "tscr"),
                             d0 = 0.5,
                             distribution.t0 = "normal",
                             initpar = NULL,
                             fixpar = NULL,
                             integrate.r = "kronrod31",
                             n_cores = 1){
  #-----------------------------------------------------------------------------
  # Set up data and model definition
  #-----------------------------------------------------------------------------
  base <- DataSetUp(data, variables = variables, correction = FALSE, history = TRUE)

  npar <- 7
  if(!is.null(initpar) & length(initpar) != npar){stop("Wrong number of parameters initialized")}
  if(is.null(initpar) & !is.null(fixpar)){stop("If some parameters are fixed please initialize the parameters using initpar.")}
  if(!is.null(fixpar) & length(initpar)!=length(fixpar)){stop("The fixpar and initpar arguments must be of the same length")}
  if(is.null(initpar)){initpar <- c(-0.841,-0.25,-8.75,-4.75,0.56,2.000,2.000)}
  if(is.null(fixpar)){fixpar <- rep(F, length(initpar))}
  #-----------------------------------------------------------------------------
  # Model Likelihood
  #-----------------------------------------------------------------------------
  nb_iter <- NULL
  nb_iter <- 0
  cat("Estimating...\n")
  start <- Sys.time()
  model <- list(par = initpar, fix = fixpar, convergence = FALSE)
  model$distribution <- list(onset = distribution.t0)
  v0 <- Volume(d0)
  #-----------------------------------------------------------------------------
  LogLikelihood <- function(par,base,distrib,v0,d0,integrate.r,n_cores = 1,initial = FALSE){
    #---------------------------------------------------------------------------
    # Parameters
    #---------------------------------------------------------------------------
    p <- ParSetUp(base, par, distrib, d0)
    #---------------------------------------------------------------------------
    # Integrand functions
    #---------------------------------------------------------------------------
    integrand_scr <- function(r,v,tdet,tscr,scr,p,v0,d0,distrib){
      if(length(scr) > 0){
        pscr <- pmin(pmax(Sensitivity(Diameter(v * exp((tscr - tdet) / r)), p[4], p[5], d0), 1e-16), 1 - 1e-16)
        phist <- exp(sum(ifelse(scr == 1, log(pscr), 0) + ifelse(scr == 0, log(1 - pscr), 0)))
      }else{phist <- 1}
      return(
        phist *
          exp(-p[3]*r*(v-v0)) * r/v *
          pdf_r(r,p[1],p[2]) *
          pdf_t0(tdet - r*log(v/v0),p[6],p[7],distrib$onset)
      )
    }
    integrand_scr <- Vectorize(integrand_scr,c('r'))

    integrand_sym <- function(r,v,tdet,tscr,scr,p,v0,d0,distrib){
      if(length(scr) > 0){
        pscr <- pmin(pmax(Sensitivity(Diameter(v * exp((tscr - tdet) / r)), p[4], p[5], d0), 1e-16), 1 - 1e-16)
        phist <- exp(sum(ifelse(scr == 1, log(pscr), 0) + ifelse(scr == 0, log(1 - pscr), 0)))
      }else{phist <- 1}
      return(
        phist *
          dexp(v-v0, p[3]*r) *
          pdf_r(r,p[1],p[2]) *
          pdf_t0(tdet - r*log(v/v0),p[6],p[7],distrib$onset)
      )
    }
    integrand_sym <- Vectorize(integrand_sym,c('r'))
    #---------------------------------------------------------------------------
    # Likelihood calculation
    #---------------------------------------------------------------------------
    n <- nrow(base)

    lld <- foreach::foreach(i = itertools::isplitVector(1:n, chunks = n_cores),
                            .combine = c
    ) %dopar% {
      vapply(i, function(i){
        p_i <- unlist(sapply(p, function(x){if(length(x) > 1){return(x[i])}else{return(x)}}))
        det <- base$det[i]
        tdet <- base$tdet[i]
        scr <- base$scr[[i]]
        tscr <- base$tscr[[i]]

        if(det %in% c(1,2)){
          ddet <- base$ddet[i]
          vdet <- Volume(ddet)

          r_low <- 0
          r_up <- tdet/log(vdet/v0)
          weights <- 0.5 * (r_up - r_low) * kronrod_w
          nodes <- 0.5 * (r_up - r_low) * kronrod_x + 0.5 * (r_up + r_low)
          if(det == 1){
            if(integrate.r == "kronrod31"){
              res <- sum(weights * integrand_sym(nodes,vdet,tdet,tscr,scr,p_i,v0,d0,distrib))
            }else if(integrate.r == "adaptative"){
              res <- integrate(integrand_sym, r_low, r_up, v=vdet, tdet=tdet, tscr=tscr, scr=scr, p=p_i, v0=v0, d0=d0, distrib=distrib, stop.on.error = FALSE)$value
            }
          }else if(det == 2){
            if(integrate.r == "kronrod31"){
              res <- sum(weights * integrand_scr(nodes,vdet,tdet,tscr,scr,p_i,v0,d0,distrib))
            }else if(integrate.r == "adaptative"){
              res <- integrate(integrand_scr, r_low, r_up, v=vdet, tdet=tdet, tscr=tscr, scr=scr, p=p_i, v0=v0, d0=d0, distrib=distrib, stop.on.error = FALSE)$value
            }
          }
        }
        res
      }, FUN.VALUE = numeric(1))
    }
    if(initial == TRUE){return(log(lld))}else{return(sum(log(lld)))}
  }
  #-----------------------------------------------------------------------------
  # The function to optimize to estimate model parameters
  fct_optim <- function(p){
    nb_iter <<- nb_iter + 1
    model$par[!model$fix] <- p # initpar value change to p if fixpar is F
    LogL <- LogLikelihood(model$par, base, model$distribution, v0, d0, integrate.r, n_cores)
    LogL <- ifelse(is.nan(LogL)|is.na(LogL)|!is.finite(LogL), -1e16, LogL)
    return(LogL)
  }
  #-----------------------------------------------------------------------------
  i <- NULL # dumb R CMD check
  cl = parallel::makeCluster(n_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  # Check of the first likelihood evaluation
  check <- LogLikelihood(model$par, base, model$distribution, v0, d0, integrate.r, n_cores, initial = TRUE)
  if(sum(is.nan(check))>0){stop("Likelihood cannot be assessed for patients ", paste(base$id[which(is.nan(check))],"")," please try other initial parameters values")}
  if(sum(!is.finite(check))>0){stop("Infinite likelihood due to patients ", paste(base$id[which(!is.finite(check))],"")," please try other initial parameters values")}
  # Optimization
  est <- marqLevAlg::mla(b = model$par[!model$fix], fn = fct_optim, minimize = FALSE, epsa = 1e-3, epsb = 1e-3, epsd = 1e-3, maxiter = 50, print.info = FALSE, blinding = TRUE) #Optimization call method
  parallel::stopCluster(cl)
  #-----------------------------------------------------------------------------
  # Set up of the output
  #-----------------------------------------------------------------------------
  model$par[!model$fix] <- est$b
  model$lld <- est$fn.value
  hessian <- diag(length(model$par[!model$fix]))
  hessian[upper.tri(hessian, diag=TRUE)] <- est$v
  model$vcov <- hessian + t(hessian) - diag(diag(hessian))
  if(est$istop == 1){
    model$convergence <- TRUE
    model$iterations <- est$ni
    model$time <- round(Sys.time()-start,3)
    model$aic <- 2*length(model$par[!model$fix]) - 2*model$lld
  }
  model$growth <- c("-log(tau1)","log(tau1/tau2)")
  model$symptoms <- c("log(eta)")
  model$sensitivity <- c("beta0","beta1")
  if(distribution.t0 == "weibull"){model$onset <- c("log(k)","log(lambda)")
  }else if(distribution.t0 == "gompertz"){model$onset <- c("log(a)","log(b)")
  }else if(distribution.t0 == "gamma"){model$onset <- c("-log(a)","log(a/b)")
  }else if(distribution.t0 == "log-normal"){model$onset <- c("meanlog","log(sdlog)")
  }else if(distribution.t0 == "normal"){model$onset <- c("log(mu)", "log(sigma)")}
  model$diameter0 <- d0
  model$call <- match.call()

  class(model) <- "TumorGrowthModel"
  return(model)
}














# TGM <- TumorGrowthModel(data_simu[1:20,],
#                         variables=list(id = "id",
#                                        det = "det",
#                                        tdet = "tdet",
#                                        ddet = "ddet",
#                                        scr = "scr",
#                                        tscr = "tscr"),
#                         d0 = 0.5,
#                         distribution.t0 = "normal",
#                         initpar = NULL,
#                         fixpar = NULL,
#                         integrate.r = "kronrod31",
#                         n_cores = 10)
