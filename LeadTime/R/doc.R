#' LeadTime: R package for lead time bias estimation and correction in cancer screening studies.
#'
#' @description
#' LeadTime is a package designed to analyze cancer screening data and estimate
#' lead time bias. This bias is due to the difference in temporality between a
#' diagnosis by screening and a diagnosis by symptoms. Consequently, the observed
#' post-diagnosis survival is overestimated for patients detected by screening
#' compared to those detected by symptoms. The cause is that the survival times
#' compared have different starting points in the natural history of the disease
#' because screening always takes place before the appearance of symptoms. The
#' longer survival observed in patients is partly due to the benefit of early
#' detection of cancer but is also artificially extended by lead time bias.
#'
#' @author Marius Robert
#'
#' @section LeadTime functions:
#' \code{\link{LeadTimeCor}}
#' \code{\link{TumorGrowthModel}}
#'
#' @keywords internal
"_PACKAGE"
NULL
