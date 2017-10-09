#########################
### global constants ####
#########################

#' bigleaf constants
#' 
#' @description Constants used in the bigleaf package
#'
#'
#' @export
bigleaf.constants <- function(){
  
  list(
    cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
    Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
    Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull_1988 p.641)
    Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken p. 245)
    Md         = 0.0289645,       # molar mass of dry air [kg mol-1]
    Mw         = 0.0180153,       # molar mass of water vapor [kg mol-1] 
    eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (-) (=Mw/Md)
    Kelvin     = 273.15,          # conversion degree Celsius to Kelvin
    g          = 9.81,            # gravitational acceleration (m s-2)
    pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
    Tair0      = 273.15,          # reference air temperature (K)
    k          = 0.41,            # von Karman constant (-)
    Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
    Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
    sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
    DwDc       = 1.6,             # Ratio of the molecular diffusivities for water vapor and CO2 (-)
    Rbwc       = 1.37             # Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
  )
  
}