#' Spherical tumor volume to diameter conversion function
#'
#' @param v numerical term for the volume of a spherical tumor
#'
#' @return the diameter of a spherical tumour of volume v
#' @export
Diameter <- function(v){
  return( ((6/pi)*v)^(1/3) )
}



#' Spherical tumor diameter to volume conversion function
#'
#' @param d numerical term for the diameter of a spherical tumor
#'
#' @return the volume of tumour of diameter d
#' @export
Volume <- function(d){
  return( (pi/6)*(d^3))
}
