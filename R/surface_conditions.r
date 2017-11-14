#############################
#### Surface conditions  ####
#############################

#' Big-leaf surface conditions
#' 
#' @description Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' @param data       Data.frame or matrix containing all required input variables
#' @param Tair       Air temperature (deg C)
#' @param pressure   Atmospheric pressure (kPa)
#' @param H          Sensible heat flux (W m-2)
#' @param LE         Latent heat flux (W m-2)
#' @param VPD        Vapor pressure deficit (kPa)
#' @param Ga         Aerodynamic conductance for heat and water vapor (m s-1)
#' @param calc.Csurf Calculate surface CO2 concentration?
#' @param Ca         Atmospheric CO2 concentration (mol mol-1). Required if calc.Ca = TRUE
#' @param NEE        Net ecosystem exchange (umol m-2 s-1). Required if calc.Ca = TRUE
#' @param Ga_CO2     Aerodynamic conductance for CO2 (m s-1). Required if calc.Ca = TRUE          
#' @param constants  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr 
#'                   eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#' 
#' @details Canopy surface temperature and humidity are calculated by inverting bulk transfer equations of
#'          sensible and latent heat, respectively. 'Canopy surface' in this case refers to 
#'          the surface of the big-leaf (i.e. at height d + z0h; the apparent sink of sensible heat and water vapor).
#'          Aerodynamic canopy surface temperature is given by:
#'          
#'          \deqn{Tsurf = Tair + H / (\rho * cp * Ga)}
#'          
#'          Vapor pressure at the canopy surface is:
#'          
#'          \deqn{esurf = e + (LE * \gamma)/(Ga * \rho * cp)}
#'          
#'          Vapor pressure deficit (VPD) at the canopy surface is calculated as:
#'          
#'          \deqn{VPD_surf = Esat_surf - esurf}
#'          
#'          Note that Ga is assumed to be equal for water vapor and sensible heat.
#'          Ga is further assumed to be the inverse of the sum of the turbulent part
#'          and the canopy boundary layer conductance (1/Ga = 1/Ga_m + 1/Gb; 
#'          see \code{\link{aerodynamic.conductance}}). If Ga is replaced by Ga_m (i.e.
#'          only the turbulent conductance part), the results of the functions represent
#'          conditions outside the canopy boundary layer, i.e. in the canopy airspace.
#' 
#' @return a data.frame with the following columns:
#'         \item{Tsurf}{Surface temperature (deg C)} \cr
#'         \item{esat_surf}{Saturation vapor pressure at the surface (kPa)} \cr
#'         \item{esurf}{vapor pressure at the surface (kPa)} \cr
#'         \item{VPD_surf}{vapor pressure deficit at the surface (kPa)} \cr
#'         \item{qsurf}{specific humidity at the surface (kg kg-1)} \cr
#'         \item{rH_surf}{relative humidity at the surface (-)} \cr
#'         \item{Ca_surf}{CO2 concentration at the surface (umol mol-1)}             
#'         
#' @examples
#' # calculate surface temperature, water vapor, VPD etc. at the surface
#' # for a given temperature and turbulent fluxes, and under different 
#' # aerodynamic conductance.
#' surface.conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1)) 
#'          
#' # now calculate also surface CO2 concentration
#' surface.conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1),
#'                    Ca=400,Ga_CO2=c(0.02,0.05,0.1),NEE=-20,calc.Csurf=TRUE)
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#' @export 
surface.conditions <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",
                               VPD="VPD",Ga="Ga",calc.Csurf=FALSE,Ca="Ca",Ga_CO2="Ga_CO2",
                               NEE="NEE",constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,LE,H,VPD,Ga))
  
  if (calc.Csurf){
    check.input(data,list(Ca,NEE,Ga_CO2))
  }
  
  rho   <- air.density(Tair,pressure)
  gamma <- psychrometric.constant(Tair,pressure,constants)
  
  # 1) Temperature
  Tsurf <- Tair + H / (rho * constants$cp * Ga)
  
  # 2) Humidity
  esat      <- Esat(Tair)[,"Esat"]
  e         <- esat - VPD
  esat_surf <- Esat(Tsurf)[,"Esat"]
  esurf     <- e + (LE * gamma)/(Ga * rho * constants$cp)
  VPD_surf  <- pmax(esat_surf - esurf,0)
  qsurf     <- VPD.to.q(VPD_surf,Tsurf,pressure,constants)
  rH_surf   <- VPD.to.rH(VPD_surf,Tsurf)
  
  # 3) CO2 concentration
  Ca_surf <- as.numeric(rep(NA,length(Tsurf)))
  if (calc.Csurf){
    Ca_surf <- Ca.surface(Ca,NEE,Ga_CO2,Tair,pressure)
  }
  
  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
}



#' CO2 concentration at the canopy surface
#'
#' @description the CO2 concentration at the canopy surface derived from net ecosystem
#'              CO2 exchange and measured atmospheric CO2 concentration.
#'              
#' @param Ca       Atmospheric CO2 concentration (umol mol-1)
#' @param NEE      Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param Ga_CO2   Aerodynamic conductance for CO2 (m s-1)
#' @param Tair     Air temperature (degC)
#' @param pressure Atmospheric pressure (kPa)
#' 
#' @details CO2 concentration at the canopy surface is calculated as:
#' 
#'        \deqn{Ca_surf = Ca + NEE / Ga_CO2}
#'        
#'        Note that this equation can be used for any gas measured (with NEE
#'        replaced by the net exchange of the respective gas and Ga_CO2 by the Ga of 
#'        that gas).
#' 
#' @note the following sign convention is employed: negative values of NEE denote
#'       net CO2 uptake by the ecosystem.
#' 
#' @return \item{Ca_surf -}{CO2 concentration at the canopy surface (umol mol-1)}
#' 
Ca.surface <- function(Ca,NEE,Ga_CO2,Tair,pressure){
  
  Ga_CO2 <- ms.to.mol(Ga_CO2,Tair,pressure)
  
  Ca_surf <- Ca + NEE/Ga_CO2
  
  return(Ca_surf)
  
}



#' Radiometric surface temperature
#' 
#' @description Radiometric surface temperature from longwave upward radiation
#'              measurements.
#'              
#' @param longwave.up longwave upward radiation (W m-2)
#' @param emissivity  infrared emissivity of the surface (-)
#' @param constants   sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                    Kelvin - conversion degree Celsius to Kelvin 
#' 
#' @details Radiometric surface temperature (Trad) is calculated as:
#' 
#'          \deqn{Trad = LW_up / (\sigma \epsilon)^(1/4)}   
#' 
#' @return a data.frame with the following columns:
#'         \item{Trad_K}{Radiometric surface temperature (K)} \cr
#'         \item{Trad_degC}{Radiometric surface temperature (degC)} 
#' 
#' @examples 
#' # determine radiative temperature of an object that has an emissivity of 0.98 
#' # and emits longwave radiation of 400Wm-2  
#' Trad.surface(400,0.98)
#' 
#' @export
Trad.surface<- function(longwave.up,emissivity,constants=bigleaf.constants()){
  
  Trad.K    <- (longwave.up / (constants$sigma * emissivity))^(1/4)
  Trad.degC <- Trad.K - constants$Kelvin
  
  return(data.frame(Trad.K,Trad.degC))
}