#' Lead time bias estimation and correction using a tumor growth model.
#'
#' @param data A data frame containing the observations and variables necessary
#' for lead time bias estimation.
#' @param variables A list containing the names of variables in the inputed data frame \code{data}.
#' \describe{
#'  \item{id}{A numeric id for each subject.}
#'  \item{det}{The cancer detection mode: 1 through symptoms, 2 through screening.}
#'  \item{tdet}{The age at cancer detection.}
#'  \item{ddet}{The tumor diameter (mm) at cancer detection (not mandatory if \code{method = "Duffy"}).}
#'  \item{scr}{A list of screening results, 0 for negative, 1 for positive.}
#'  \item{tscr}{A list for each patient of their age during their previous negative screening
#'  (not mandatory  if there is \code{par} or \code{model} in input).}
#'  \item{last}{The patient's vital status on the date of latest news: 0 alive, 1 dead.}
#'  \item{tlast}{The patient's age on the date of latest news.}
#' }
#' @param par An optional vector for lead time estimation parameters. Useful for
#' bias estimation with parameters already estimated or from the literature.
#' @param model An optional model object. It can be of class \code{TriStateModel} if
#' \code{method = "Duffy"} or \code{TumorGrowthModel} if \code{method = "Abrahamsson"}. Useful if a model
#' has already been estimated on the population to save computational time.
#' @param overdiagnosis Boolean term. True, allows for an estimated lead time greater than the observed survival time for a deceased patient.
#' False, truncates the lead time distribution so that its expected value does not exceed the observed survival time.
#' @param d0 The smallest diameter of a tumor detectable by screening to pass in a tumor growth model estimation.
#' @param distribution.t0 A str term to define the distribution of age at onset.
#' @param initpar An optional vector to initialize model parameters, if a model has to be estimated.
#' @param fixpar An optional boolean vector with TRUE to fix model parameters
#' @param n_cores The number cores for parallel computation. Default to 1 (sequential mode),
#' parallelization is useful for saving time on the often long estimation of continuous tumor growth models.
#'
#' @export
LeadTimeCor <- function(data,
                        variables=list(id = "id",
                                       det = "det",
                                       tdet = "tdet",
                                       ddet = "ddet",
                                       scr = "scr",
                                       tscr = "tscr",
                                       last = "last",
                                       tlast = "tlast"),
                        par = NULL,
                        model = NULL,
                        overdiagnosis = FALSE,
                        d0 = 0.5,
                        distribution.t0 = "normal",
                        initpar = NULL,
                        fixpar = NULL,
                        n_cores = 1){
  #-----------------------------------------------------------------------------
  # Set up data and parameters
  #-----------------------------------------------------------------------------
  if(!is.null(par) & !is.null(model)){stop("The function takes input either parameters or a model but not both")}

  base <- DataSetUp(data, variables = variables, correction = TRUE)

  if(!is.null(model) & !is.null(par)){warning("par values ignored, priority on model estimations")}
  if(!is.null(model) & as.character(class(model)) == "TumorGrowthModel"){
    par <- model$par
  }
  if(!is.null(model) & as.character(class(model)) != "TumorGrowthModel"){stop("The model is not tumor growth type.")}

  if(!is.null(par) & length(par) != 7){stop("Incorrect number of parameters")}
  if(is.null(par)){
    cat("\nRunning continuous tumor growth model to estimate parameters...")
    model <- TumorGrowthModel(data,
                              variables=list(id = variables$id,
                                             det = variables$det,
                                             tdet = variables$tdet,
                                             ddet = variables$ddet,
                                             scr = variables$scr,
                                             tscr = variables$tscr),
                              d0 = d0,
                              distribution.t0 = "normal",
                              initpar=initpar,
                              fixpar=fixpar,
                              n_cores = n_cores)
    par <- model$par
  }

  n <- nrow(base)
  v0 <- Volume(d0)
  start <- Sys.time()
  #-----------------------------------------------------------------------------
  # Parameters
  #-----------------------------------------------------------------------------
  if(!is.null(model) & as.character(class(model)) == "TumorGrowthModel"){
    distrib <- model$distribution
  }else{
    distrib <- list(onset = distribution.t0)
  }
  base_scr <- subset(base, det == 2)
  nscr <- nrow(base_scr)
  p <- ParSetUp(base_scr, par, distrib)
  #-----------------------------------------------------------------------------
  # Integrand function
  #-----------------------------------------------------------------------------
  integrand_lt <- function(r,l,vdet,tdet,tscr,scr,p,v0,d0,distrib){
    if(length(scr) > 0){
      pscr <- pmin(pmax(Sensitivity(Diameter(vdet * exp((tscr - tdet) / r)), p[4], p[5], d0), 1e-16), 1 - 1e-16)
      phist <- exp(sum(ifelse(scr == 1, log(pscr), 0) + ifelse(scr == 0, log(1 - pscr), 0)))
    }else{phist <- 1}
    return(
      p[3]*vdet*exp(-p[3]*r*(vdet*exp(l/r) - vdet)+l/r) *
        phist *
        exp(-p[3]*r*(vdet-v0)) * r/vdet *
        pdf_r(r,p[1],p[2]) *
        pdf_t0(tdet - r*log(vdet/v0),p[6],p[7],distrib$onset)
    )
  }
  integrand_lt <- Vectorize(integrand_lt, c('r'))

  density_lt <- function(l,vdet,tdet,tscr,scr,p,v0,d0,distrib){
    r_low <- 0
    r_up <- tdet/log(vdet/v0)
    tol <- ifelse(vdet > 100,0.001,1/10^(3+max(which(vdet<=c(100,50,20,10,5,2,1,0.5,0.2,0.1)))))
    return(
      integrate(integrand_lt,lower=r_low,upper=r_up,l=l,vdet=vdet,tdet=tdet,tscr=tscr,scr=scr,p=p,v0=v0,d0=d0,distrib=distrib,subdivisions = 100,rel.tol = tol,stop.on.error = FALSE)$value
    )
  }
  density_lt <- Vectorize(density_lt, c('l'))
  #-----------------------------------------------------------------------------
  # Lead time density expected value
  #-----------------------------------------------------------------------------
  expected_lt <- function(det,tdet,ddet,vdet,scr,tscr,last,tlast,overdiagnosis,p,v0,d0,distrib){

    if(last == 0 | overdiagnosis == TRUE){
      try_l <- c(10,20,50,100,200,500,1000)
      l_low <- 0
      l_up <- try_l[min(which(density_lt(try_l,vdet,tdet,tscr,scr,p,v0,d0,distrib)<1e-16))]
      weights <- 0.5 * (l_up - l_low) * kronrod_w
      nodes <- 0.5 * (l_up - l_low) * kronrod_x + 0.5 * (l_up + l_low)
      numer <- sum(weights * nodes*density_lt(nodes,vdet,tdet,tscr,scr,p,v0,d0,distrib))
      denom <- sum(weights * density_lt(nodes,vdet,tdet,tscr,scr,p,v0,d0,distrib))
    }else if(last == 1 & overdiagnosis == FALSE){
      l_low <- 0
      l_up <- tlast - tdet
      weights <- 0.5 * (l_up - l_low) * kronrod_w
      nodes <- 0.5 * (l_up - l_low) * kronrod_x + 0.5 * (l_up + l_low)
      numer <- sum(weights * nodes*density_lt(nodes,vdet,tdet,tscr,scr,p,v0,d0,distrib))
      denom <- sum(weights * density_lt(nodes,vdet,tdet,tscr,scr,p,v0,d0,distrib))
    }
    return(numer/denom)
  }
  #-----------------------------------------------------------------------------
  # Parallelization, set up of the loop
  #-----------------------------------------------------------------------------
  cat("\nCorrection estimation...")
  start <- Sys.time()
  i <- NULL # dumb R CMD check
  cl = parallel::makeCluster(n_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  '%dopar%' <- foreach::'%dopar%'
  #-----------------------------------------------------------------------------
  # Correction estimation
  #-----------------------------------------------------------------------------
  corr <- foreach::foreach(i = 1:nscr, .combine = c
  ) %dopar% {
    vapply(i, function(i){
      p_i <- unlist(sapply(p, function(x){if(length(x) > 1){return(x[i])}else{return(x)}}))
      tdet <- base_scr$tdet[i]
      ddet <- base_scr$ddet[i]
      vdet <- Volume(ddet)
      scr <- base_scr$scr[[i]]
      tscr <- base_scr$tscr[[i]]
      last <- base_scr$last[i]
      tlast <- base_scr$tlast[i]

      expected_lt(det,tdet,ddet,vdet,scr,tscr,last,tlast,overdiagnosis,p_i,v0,d0,distrib)

    }, FUN.VALUE = numeric(1))
  }
  #-----------------------------------------------------------------------------
  parallel::stopCluster(cl)

  cat("\nCompleted, lead time estimated for",sum(ifelse(corr[!is.na(corr)&!is.nan(corr)]>0,1,0)), "individuals.")
  if(sum(is.na(corr)|is.nan(corr)) > 0){
    cat("Warning: estimation failed for",sum(is.na(corr)|is.nan(corr)), "individuals. \n")
  }
  cat("\nThe estimated correction is greater than the observed survival for", length(corr[corr > (base_scr$tlast - base_scr$tdet)]), "individuals.")
  if(overdiagnosis == TRUE){
    cat("\nOverdiagnosis estimated for",length(corr[base_scr$last == 1 & corr > (base_scr$tlast - base_scr$tdet)]), "individuals.")
  }
  #-----------------------------------------------------------------------------
  # Set up of the output
  #-----------------------------------------------------------------------------
  output <- list()
  output$model <- NULL
  if(!is.null(model)){output$model <- model}else{
    output$parameters <- par
  }
  output$id <- base$id
  output$detection <- base$det
  output$observed <- base$tlast - base$tdet
  output$leadtime <- ifelse(output$detection == 1, 0, NA)
  output$leadtime[output$detection == 2] <- corr
  output$corrected <- ifelse(output$observed - output$leadtime > 0, output$observed - output$leadtime, 0)
  output$death <- base$last
  output$overdetected <- ifelse(output$observed <= output$leadtime, 1, 0)
  if(overdiagnosis == TRUE){
    output$overdiagnosis <- ifelse(output$death == 1 & output$observed <= output$leadtime, 1, 0)
  }

  df <- data.frame(id = output$id,
                   detection = output$detection,
                   observed = output$observed,
                   leadtime = output$leadtime,
                   death = output$death)

  # survfit object
  survival <- list()
  if(length(df$detection[df$detection == 1]) > 1){
    sym <- subset(df, df$detection==1)
    survival$sym <- survival::survfit(survival::Surv(sym$observed,sym$death)~1)
  }
  scr <- subset(df, df$detection==2)
  survival$scr <- survival::survfit(survival::Surv(scr$observed,scr$death)~1)
  survival$cor <- survival::survfit(survival::Surv(ifelse(scr$observed - scr$leadtime > 0,scr$observed - scr$leadtime,0), scr$death)~1)
  output$survival <- survival

  output$time <- round(Sys.time()-start,3)
  class(output) <- "LeadTime"
  return(output)
}
