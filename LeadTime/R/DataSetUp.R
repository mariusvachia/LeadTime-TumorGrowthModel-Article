DataSetUp <- function(data,
                      variables=list(id = "id",
                                     det = "det",
                                     tdet = "tdet",
                                     ddet = "ddet",
                                     scr = "scr",
                                     tscr = "tscr",
                                     last = "last",
                                     tlast = "tlast"),
                      correction = FALSE,
                      history = TRUE){
  cat("Data checking and preparation... \n")
  n <- nrow(data)
  base <- data.frame(matrix(nrow = n,ncol = 0))
  has_missing <- function(x) {
    is.null(x) || any(is.na(x)) || (is.numeric(x) && any(is.nan(x)))
  }
  #-----------------------------------------------------------------------------
  # id, det, tdet
  if(is.null(variables$id) & is.null(data$id)){base$id <- c(1:n)
  }else{base$id <- data[,ifelse(is.null(variables$id),"id",variables$id)]}
  if(is.null(variables$det) & is.null(data$det)){stop("The variable det for the detection mode is required")}
  if(is.null(variables$tdet) & is.null(data$tdet)){stop("The variable tdet for the time at cancer detection is required")}
  base$det <- data[,ifelse(is.null(variables$det),"det",variables$det)]
  base$tdet <-  data[,ifelse(is.null(variables$tdet),"tdet",variables$tdet)]
  if(has_missing(base$id)){warning("Missing values for the id")}
  if(length(unique(base$id)) != n){warning("Multiple values for the id variable")}
  if(has_missing(base$det)){stop("Missing value for detection mode")}
  if(!(class(base$tdet) %in% c("numeric","integer"))){base$tdet <- as.numeric(base$tdet)}
  if(has_missing(base$tdet)){stop("tdet variable for time at cancer detection has to be numeric without missing values")}
  #-----------------------------------------------------------------------------
  # ddet, vdet
  if(is.null(variables$ddet) & is.null(data$ddet)){stop("the variable ddet for the tumor diameter at cancer detection is required")}
  base$ddet <- data[,ifelse(is.null(variables$ddet),"ddet",variables$ddet)]
  if(!(class(base$ddet) %in% c("numeric","integer"))){base$ddet <- as.numeric(base$ddet)}
  if(has_missing(base$ddet)){stop("ddet variable for tumor diameter at detection has to be numeric without missing values")}
  base$vdet <- Volume(base$ddet)
  #-----------------------------------------------------------------------------
  # last, tlast
  if(correction == TRUE){
    if(is.null(variables$last) & is.null(data$last)){stop("the variable last for the vital status at latest news is required")}
    if(is.null(variables$tlast) & is.null(data$tlast)){stop("the variable tlast for the time at latest news is required")}
    base$last <- data[,ifelse(is.null(variables$last),"last",variables$last)]
    base$tlast <- data[,ifelse(is.null(variables$tlast),"tlast",variables$tlast)]
    if(has_missing(base$last)){stop("Missing values for the vital status at latest news, variable last ")}
    if(!(class(base$tlast) %in% c("numeric","integer"))){base$tlast <- as.numeric(base$tlast)}
    if(has_missing(base$tlast)){stop("tlast variable for time at latest news has to be numeric without missing values")}
  }
  #-----------------------------------------------------------------------------
  # scr, tscr
  if(history == TRUE){
    if(is.null(variables$tscr) & is.null(data$tscr)){stop("tscr the time of the previous screening is required")}
    tscr <- data[,ifelse(is.null(variables$tscr),"tscr",variables$tscr)]
    if(has_missing(tscr)){for(i in which(is.na(tscr))){tscr[[i]] <- numeric(0)}}
    if(!is.null(variables$scr) | ("scr" %in% names(data))){
      scr <- data[,ifelse(is.null(variables$scr),"scr",variables$scr)]
      if(has_missing(scr)){for(i in which(is.na(scr))){scr[[i]] <- numeric(0)}}
    }else{
      scr <- list()
    }
  }
  #-----------------------------------------------------------------------------
  # Checking the mammography history
  det_unique <- unique(as.factor(base$det)) # number of unique values of the variable det
  if(history == TRUE){
    scr_state <- 0
    if(length(det_unique) <= 3 & length(scr) == n){
      scr_n <- numeric(n);scr_pos <- numeric(n);scr_last <- numeric(n)
      for(i in 1:n){
        if(has_missing(tscr[[i]])){stop("Missing values for tscr")}
        if(any(!(scr[[i]] %in% c(0,1)))){stop("scr variable for screening results only takes 0 or 1 values")}
        ord <- order(tscr[[i]]); tscr[[i]] <- tscr[[i]][ord]; scr[[i]] <- scr[[i]][ord]
        scr_n[i] <- length(scr[[i]]);scr_pos[i] <- sum(scr[[i]] == 1);scr_last[i] <- ifelse(length(scr[[i]]) > 0,tail(scr[[i]],1),0)
        if(sum(tscr[[i]] > base$tdet[i])){stop("screening after cancer detection are not allowed")}}
      scr_tab <- table(scr_last[scr_n > 0],base$det[scr_n > 0])
      if(length(scr_tab) > 2*ncol(scr_tab)){stop("The supported values for the screening results are 0 for negative, 1 for positive")
      }else if(sum(scr_pos) == 0){scr_state <- 1
      }else if(length(scr_tab) == 2*sum(scr_tab>0) & sum(scr_tab[2,]>0) == 1){scr_state <- 2
      }else{scr_state <- 3; warning("The number of asymptomatic cancer detections and the number of positive screenings are different")}
    }else if(length(det_unique) <= 3 & length(scr) == 0 & length(tscr) == n){
      scr_state <- 1
      message("Assuming all the screening times concern negative screening")
    }else{stop("Issue with detection mode or/and mammography history")}
  }
  #-----------------------------------------------------------------------------
  # Determining detection mode
  det_scr <- NA
  if(scr_state == 2){det_scr <- names(which(scr_tab[2,]>0))}

  if(length(det_unique) == 1){
    if(correction == TRUE){
      message("Estimates are only made on screen-detected patients")
      base$det <- rep(2,n)
    }else if(correction == FALSE){
      stop("The variable det must have values for symptomatic and screen-detected patients for model estimation")
    }
  }else if(length(det_unique) == 2){
    if(!is.na(det_scr)){
      if(det_scr != 2){
        message("Assuming det = ", det_scr, "is for screen-detected patients")
        base$det <- ifelse(as.factor(base$det) == det_scr,2,1)
      }
    }else if(sum(as.factor(unique(sort(base$det))) == c("0","1")) == 2){
      base$det <- base$det + 1
      message(paste("Assuming det =", 1," is for screen-detected patients"))
    }else if(sum(as.factor(unique(sort(base$det))) == c("1","2")) != 2){stop("Not able to determine who are the screen-detected patients")}
  }else if(length(det_unique) == 3){
    stop("The model does not support censored individual data for cancer detection")
    if(!is.na(det_scr)){
      if(det_scr != 2){stop("det value must be 0 for censored cases, 1 for symptomatic and 2 for screen-detected")}
    }
    if(sum(as.factor(unique(sort(base$det))) == c("0","1","2")) != 3){stop("det value must be 0 for censored cases, 1 for symptomatic and 2 for screen-detected")}
  }else{stop("Too much unique values for variable det")}
  #-----------------------------------------------------------------------------
  # Assembling full history
  if(history == TRUE & length(scr) == 0 & scr_state == 1){
    for(i in 1:n){
      ord <- order(tscr[[i]]); tscr[[i]] <- tscr[[i]][ord]
      if(has_missing(tscr[[i]])){stop("Missing values for tscr")}
      scr[[i]] <- numeric(0)
      if(length(tscr[[i]]) > 0 & base$det[i] %in% c(0,1)){
        scr[[i]] <- rep(0,length(tscr[[i]]))
      }else if(length(tscr[[i]]) > 0 & base$det[i] == 2){
        scr[[i]] <- c(rep(0,length(tscr[[i]])),1)
        tscr[[i]] <- c(tscr[[i]],base$tdet[i])
      }else if(base$det[i] == 2){
        scr[[i]] <- 1
        tscr[[i]] <- base$tdet[i]
      }
    }
  }
  if(history == TRUE){
    base$scr <- scr
    base$tscr <- tscr
  }

  return(base)
}
