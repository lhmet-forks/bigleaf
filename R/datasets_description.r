#' Eddy covariance data of DE-Tha (Tharandt)
#' 
#' @description Halfhourly eddy covariance Data of the site DE-Tha,
#'              a spruce forest in Eastern Germany 
#'              (\url{http://www.fluxdata.org:8080/sitepages/siteInfo.aspx?DE-Tha}). 
#'              Data are from June 2014.
#'              
#' @format A data frame with 1440 observations and 29 columns:
#'  \describe{
#'    \item{year}{year of measurement}
#'    \item{month}{month of measurement}
#'    \item{doy}{day of year}
#'    \item{hour}{hour (0 - 23.5)}
#'    \item{Tair}{Air temperature (degC)}
#'    \item{Tair_qc}{Quality control of \code{Tair}}
#'    \item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1)}
#'    \item{PPFD_qc}{Quality control of \code{PPFD}}
#'    \item{VPD}{Vapor pressure deficit (kPa)}
#'    \item{VPD_qc}{Quality control of \code{VPD}}
#'    \item{pressure}{Atmospheric pressure (kPa)}
#'    \item{precip}{precipitation (mm)}
#'    \item{precip_qc}{Quality control of \code{precip}}
#'    \item{ustar}{friction velocity (ms-1)}
#'    \item{wind}{horizontal wind velocity (m s-1)}
#'    \item{wind_qc}{Quality control of \code{wind}}
#'    \item{Ca}{CO2 concentration (ppm)}
#'    \item{Ca_qc}{Quality control of \code{Ca}}
#'    \item{Rn}{Net radiation (Wm-2)}
#'    \item{G}{Ground heat flux (Wm-2)}
#'    \item{G_qc}{Quality control of \code{G}}
#'    \item{LE}{Latent heat flux (Wm-2)}
#'    \item{LE_qc}{Quality control of \code{LE}}
#'    \item{H}{Sensible heat flux (Wm-2)}
#'    \item{H_qc}{Quality control of \code{H}}
#'    \item{NEE}{Net ecosystem exchange (umol m-2 s-1)}
#'    \item{NEE_qc}{Quality control of \code{NEE}}
#'    \item{GPP}{Gross primary productivity from nighttime partitioning (umol m-2 s-1)}
#'    \item{GPP_qc}{Quality control of \code{GPP}}
#'  }
#'  
#' @note Squared brackets denote the original variables as provided by the FLUXNET2015 dataset.
#'       Note that some variable units have been converted (e.g. VPD from hPa to kPa).    
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet.fluxdata.org/}  
"DE_Tha_Jun_2014"



#' Eddy covariance data of FR-Pue (Puechabon)
#' 
#' @description Halfhourly eddy covariance Data of the site FR-Pue,
#'              a Mediterranean evergreen oak forest in Southern France
#'              (\url{http://www.fluxdata.org:8080/sitepages/siteInfo.aspx?FR-Pue}).
#'              Data are from May 2012.
#'              
#' @format A data frame with 1488 observations and 27 columns:
#'  \describe{
#'    \item{year}{year of measurement}
#'    \item{month}{month of measurement}
#'    \item{doy}{day of year}
#'    \item{hour}{hour (0 - 23.5)}
#'    \item{Tair}{Air temperature (degC)}
#'    \item{Tair_qc}{Quality control of \code{Tair}}
#'    \item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1)}
#'    \item{PPFD_qc}{Quality control of \code{PPFD}}
#'    \item{VPD}{Vapor pressure deficit (kPa)}
#'    \item{VPD_qc}{Quality control of \code{VPD}}
#'    \item{pressure}{Atmospheric pressure (kPa)}
#'    \item{precip}{precipitation (mm)}
#'    \item{precip_qc}{Quality control of \code{precip}}
#'    \item{ustar}{friction velocity (ms-1)}
#'    \item{wind}{horizontal wind velocity (m s-1)}
#'    \item{wind_qc}{Quality control of \code{wind}}
#'    \item{Ca}{CO2 concentration (ppm)}
#'    \item{Ca_qc}{Quality control of \code{Ca}}
#'    \item{Rn}{Net radiation (Wm-2)}
#'    \item{LE}{Latent heat flux (Wm-2)}
#'    \item{LE_qc}{Quality control of \code{LE}}
#'    \item{H}{Sensible heat flux (Wm-2)}
#'    \item{H_qc}{Quality control of \code{H}}
#'    \item{NEE}{Net ecosystem exchange (umol m-2 s-1)}
#'    \item{NEE_qc}{Quality control of \code{NEE}}
#'    \item{GPP}{Gross primary productivity from nighttime partitioning (umol m-2 s-1)}
#'    \item{GPP_qc}{Quality control of \code{GPP}}
#'  }
#'  
#' @note Squared brackets denote the original variables as provided by the FLUXNET2015 dataset.
#'       Note that some variable units have been converted (e.g. VPD from hPa to kPa).    
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet.fluxdata.org/}  
"FR_Pue_May_2012"