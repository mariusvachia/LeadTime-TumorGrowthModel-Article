ParSetUp <- function(data,
                     p,
                     distrib,
                     d0){
  n <- nrow(data)
  res <- list(p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, p5 = NULL, p6 = NULL, p7 = NULL, p8 = NULL)
  #-----------------------------------------------------------------------------
  # Tumor growth parameters
  #-----------------------------------------------------------------------------
  res$p1 <- 1/exp(p[1])
  res$p2 <- res$p1/exp(p[2])
  #-----------------------------------------------------------------------------
  # Volume at symptoms parameters
  #-----------------------------------------------------------------------------
  res$p3 <- exp(p[3])
  #-----------------------------------------------------------------------------
  # Screening sensitivity parameters
  #-----------------------------------------------------------------------------
  res$p4 <- p[4]
  res$p5 <- p[5]
  #-----------------------------------------------------------------------------
  # Age at tumor onset parameters
  #-----------------------------------------------------------------------------
  if(distrib$onset == "weibull"){
    res$p6 <- exp(p[6])
    res$p7 <- exp(p[7])
    #-----------------------------
  }else if(distrib$onset == "gompertz"){
    res$p6 <- exp(p[6])
    res$p7 <- exp(p[7])
    #-----------------------------
  }else if(distrib$onset == "gamma"){
    res$p6 <- 1/exp(p[6])
    res$p7 <- res$p6/exp(p[7])
    #-------------------
  }else if(distrib$onset == "log-normal"){
    res$p6 <- p[6]
    res$p7 <- exp(p[7])
    #-------------------
  }else if(distrib$onset == "normal"){
    res$p6 <- exp(p[6])
    res$p7 <- exp(p[7])
  }
  return(res)
}
