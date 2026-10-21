#' Plot observed and corrected survival curves for lead time.
#'
#' @param x A LeadTime class object.
#' @param conf.int Boolean term to add confidence intervals in plot.
#' @param ylim The y limits of the plot.
#' @param xlim The x limits of the plot.
#' @param ... Argument to ignore, output uses \pkg{ggplot2}.
#'
#' @export
plot.LeadTime <- function(x, conf.int = FALSE, ylim=c(0.5,1),xlim = c(0,10),...){

  leg <- c("Screen-detected")
  pal <- c("#9c005e")
  if(!is.null(x$survival$sym)){
    leg <- c("Symptomatically detected",leg)
    pal <- c("#34415e", pal)
  }

  leg <- c(leg,"Corrected lead time screen-detected survival")
  pal <- c(pal,"#1abc9c")

  ggsurv  <- survminer::ggsurvplot(x$survival,combine = TRUE,censor = FALSE,xlab="Years from diagnosis",
                                   ylab="Survival probability",ylim=ylim, break.time.by=1,
                                   conf.int = conf.int,
                                   legend = "bottom", xlim=xlim,
                                   legend.labs = leg,
                                   palette = pal,
                                   ggtheme = ggplot2::theme_bw())# Change ggplot2 theme


  gg <- ggsurv$plot + ggplot2::theme(legend.text = ggplot2::element_text(size = 11),
                                     legend.title = ggplot2::element_blank(),
                                     legend.direction = "vertical")
  gg
  return(gg)
}
