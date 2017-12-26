#################################
### Meteorological variables ####
#################################

#' Air density
#' 
#' @description Air density of moist air from air temperature and pressure.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants Kelvin - conversion degC to Kelvin \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1)
#' 
#' @details Air density (\eqn{\rho}) is calculated as:
#' \deqn{\rho = pressure / Rd * Tair}
#' 
#' @return \item{\eqn{\rho}}{air density (kg m-3)}
#' 
#' @examples 
#' # air density at 25degC and standard pressure (101.325kPa)
#' air.density(25,101.325)
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
#' @export
air.density <- function(Tair,pressure,constants=bigleaf.constants()){
  Tair     <- Tair + constants$Kelvin
  pressure <- pressure*1000
  
  rho <- pressure / (constants$Rd * Tair) 
  
  return(rho)
}



#' Atmospheric pressure from hypsometric equation
#' 
#' @description An estimate of mean pressure at a given elevation as predicted by the
#'              hypsometric equation.
#'              
#' @param elev      Elevation a.s.l. (m)
#' @param Tair      Air temperature (deg C)
#' @param VPD       Vapor pressure deficit (kPa); optional
#' @param constants Kelvin- conversion degC to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  g -  gravitational acceleration (m s-2) \cr
#' 
#' @details Atmospheric pressure is approximated by the hypsometric equation:
#' 
#'       \deqn{pressure = pressure_0 / (exp(g * elevation / (Rd Temp)))}
#'       
#' @note The hypsometric equation gives an estimate of the standard pressure
#'       at a given altitude. 
#'       If VPD is provided, humidity correction is applied and the
#'       virtual temperature instead of air temperature is used. VPD is 
#'       internally converted to specific humidity.
#'
#' @return \item{pressure -}{Atmospheric pressure (kPa)}
#'                            
#' @references Stull B., 1988: An Introduction to Boundary Layer Meteorology.
#'             Kluwer Academic Publishers, Dordrecht, Netherlands.
#' 
#' @examples
#' # mean pressure at 500m altitude at 25 deg C and VPD of 1 kPa
#' pressure.from.elevation(500,Tair=25,VPD=1)
#' 
#' @export                           
pressure.from.elevation <- function(elev,Tair,VPD=NULL,constants=bigleaf.constants()){
  
  Tair     <- Tair + constants$Kelvin
  
  if(is.null(VPD)){
    
    pressure <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tair))
    
  } else {
    
    pressure1   <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tair))
    Tv          <- virtual.temp(Tair - constants$Kelvin,pressure1 / 1000,VPD,constants) + constants$Kelvin
    
    pressure    <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tv))
  }
  
  pressure <- pressure/1000
  return(pressure)
} 



#' Saturation vapor pressure (Esat) and slope of the Esat curve
#' 
#' @description Calculates Saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' @param Tair     Air temperature (deg C)
#' @param formula  Formula to be used. Either "Sonntag_1990" (Default) or "Alduchov_1996"
#' 
#' @details Esat (kPa) is calculated using the Magnus equation:
#' 
#'  \deqn{Esat = a * exp((b * Tair) / (c + Tair)) / 1000}
#'  
#'  where the coefficients a, b, c take different values depending on the formula used.
#'  The default values are from Sonntag 1990 (a=611.2, b=17.62, c=243.12). This version
#'  of the Magnus equation is recommended by the WMO (WMO 2008; p1.4-29).
#'  The slope of the Esat curve (\eqn{\Delta}) is calculated as the first derivative of the function:
#'  
#'  \deqn{\Delta = dEsat / dTair}
#'  
#'  which is solved using \code{\link[stats]{D}}.
#' 
#' @return A dataframe with the following columns: 
#'  \item{Esat}{Saturation vapor pressure (kPa)}
#'  \item{Delta}{Slope of the saturation vapor pressure curve (kPa K-1)}
#'    
#' @references Sonntag D. 1990: Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90 and psychrometric formulae. 
#'             Zeitschrift fuer Meteorologie 70, 340-344.
#'             
#'             Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure. Journal of Applied Meteorology, 35, 601-609
#'             
#'             World Meteorological Organization 2008: Guide to Meteorological Instruments
#'             and Methods of Observation (WMO-No.8). World Meteorological Organization,
#'             Geneva. 7th Edition.
#' 
#' @examples 
#' Esat(seq(0,45,5))[,"Esat"]  # Esat in kPa
#' Esat(seq(0,45,5))[,"Delta"] # the corresponding slope of the Esat curve in kPa K-1
#'        
#' @importFrom stats D                  
#' @export
Esat <- function(Tair,formula=c("Sonntag_1990","Alduchov_1996")){
  
  formula <- match.arg(formula)
  
  if (formula == "Sonntag_1990"){
    a <- 611.2
    b <- 17.62
    c <- 243.12
  } else if (formula == "Alduchov_1996"){
    a <- 610.94
    b <- 17.625
    c <- 243.04
  }
  
  # saturation vapor pressure
  Esat <- a * exp((b * Tair) / (c + Tair))
  Esat <- Esat/1000
  
  # slope of the saturation vapor pressure curve
  Delta <- eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  Delta <- Delta/1000
  
  
  return(data.frame(Esat,Delta))
}



#' Psychrometric constant
#' 
#' @description Calculates the psychrometric 'constant'.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' @details The psychrometric constant (\eqn{\gamma}) is given as:
#' 
#'  \deqn{\gamma = cp * pressure / (eps * \lambda)},
#'  
#'  where \eqn{\lambda} is the latent heaf of vaporization (J kg-1), 
#'  as calculated with \code{\link{latent.heat.vaporization}}.
#'  
#' @return \item{\eqn{\gamma} -}{the psychrometric constant (kPa K-1)}
#'  
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#' 
#' @examples 
#' psychrometric.constant(seq(5,45,5),100)
#' 
#' @export
psychrometric.constant <- function(Tair,pressure,constants=bigleaf.constants()){
  
  lambda <- latent.heat.vaporization(Tair)
  gamma  <- (constants$cp * pressure) / (constants$eps * lambda)
  
  return(gamma)
}



#' Latent heat of vaporization
#' 
#' @description Latent heat of vaporization as a function of air temperature.
#' 
#' @param Tair  Air temperature (deg C)
#' 
#' @details The following formula is used:
#' 
#' \deqn{\lambda = (2.501 - 0.00237*Tair)10^6}
#' 
#' @return \item{\eqn{\lambda} -}{Latent heat of vaporization (J kg-1)} 
#' 
#' @references Stull, B., 1988: An Introduction to Boundary Layer Meteorology (p.641)
#'             Kluwer Academic Publishers, Dordrecht, Netherlands
#' 
#' @examples 
#' latent.heat.vaporization(seq(5,45,5))             
#'             
#' @export
latent.heat.vaporization <- function(Tair) {
  k1 <- 2.501
  k2 <- 0.00237
  lambda <- ( k1 - k2 * Tair ) * 1e+06
  
  return(lambda)
}



#' Wet-bulb temperature
#' 
#' @description calculates the wet bulb temperature, i.e. the temperature
#'              that the air would have if it was saturated.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' @details Wet-bulb temperature (Tw) is calculated from the following expression:
#'          
#'          \deqn{e = Esat(Tw) - gamma* (Tair - Tw)}
#'          
#'          The equation is solved for Tw using \code{\link[stats]{nls}}.
#'          Actual vapor pressure e is calculated from VPD using the function \code{\link{VPD.to.e}}.
#'          
#' @return \item{Tw -}{wet-bulb temperature (degC)}      
#'              
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'             
#' @examples 
#' wetbulb.temp(Tair=c(20,25),pressure=100,VPD=c(1,1.6))             
#'        
#' @importFrom stats nls                  
#' @export
wetbulb.temp <- function(Tair,pressure,VPD,constants=bigleaf.constants()){
  
  gamma  <- psychrometric.constant(Tair,pressure)
  ea     <- VPD.to.e(VPD,Tair)
  
  wetbulb <- function(ea,gamma,Tair){
    mod <- nls(ea ~ Esat(Tw)[,"Esat"] - 0.93*gamma*(Tair - Tw),start=list(Tw=Tair),algorithm="port")
    return(summary(mod)$coef[1,1])
  }
  
  Tw <- sapply(seq_along(ea),function(i) wetbulb(ea[i],gamma[i],Tair[i]))
  
  return(Tw)
  
}



#' Dew point
#' 
#' @description calculates the dew point, the temperature to which air must be 
#'              cooled to become saturated (i.e. e = esat(Td))
#'
#' @param Tair     Air temperature (degC)
#' @param VPD      Vapor pressure deficit (kPa)
#' 
#' @details Dew point temperature (Td) is defined by:
#' 
#'         \deqn{e = Esat(Td)}
#'    
#'          which is solved for Td using \code{\link[stats]{nls}}.
#'          
#' @return \item{Td -}{dew point temperature (degC)}
#' 
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'             
#' @examples
#' dew.point(25,1.5)                
#' 
#' @importFrom stats nls 
#' @export              
dew.point <- function(Tair,VPD){
  
  ea     <- VPD.to.e(VPD,Tair)
  
  dewpoint <- function(ea,Tair){
    mod <- nls(ea ~ Esat(Td)[,"Esat"],start=list(Td=Tair),algorithm="port")
    return(summary(mod)$coef[1,1])
  }
  
  Td <- sapply(seq_along(ea),function(i) dewpoint(ea[i],Tair[i]))
  
  return(Td)
}



#' Virtual temperature
#' 
#' @description Virtual temperature, defined as the temperature at which dry air would have the same
#'              density as moist air at its actual temperature.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' @details the virtual temperature is given by:
#'  
#'  \deqn{Tv = Tair / (1 - (1 - eps) e/pressure)}
#' 
#'  where Tair is in Kelvin (converted internally). Likewise, VPD is converted 
#'  to actual vapor pressure (e in kPa) with \code{\link{VPD.to.e}} internally.
#' 
#' @return \item{Tv -}{virtual temperature (deg C)}
#' 
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'  
#' @examples 
#' virtual.temp(25,100,1.5)                        
#'               
#' @export
virtual.temp <- function(Tair,pressure,VPD,constants=bigleaf.constants()){
  
  e    <- VPD.to.e(VPD,Tair)
  Tair <- Tair + constants$Kelvin
  
  Tv <- Tair / (1 - (1 - constants$eps) * e/pressure) 
  Tv <- Tv - constants$Kelvin
  
  return(Tv)
}



#' Kinematic viscosity of air
#' 
#' @description calculates the kinematic viscosity of air.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K)
#' 
#' @details where v is the kinematic viscosity of the air (m2 s-1), 
#'          given by (Massman 1999b):
#'          
#'          \deqn{v = 1.327 * 10^-5(pressure0/pressure)(Tair/Tair0)^1.81}
#'          
#' @return \item{v -}{kinematic viscosity of air (m2 s-1)}
#' 
#' @references Massman, W.J., 1999b: Molecular diffusivities of Hg vapor in air, 
#'             O2 and N2 near STP and the kinematic viscosity and thermal diffusivity
#'             of air near STP. Atmospheric Environment 33, 453-457.      
#'             
#' @examples 
#' kinematic.viscosity(25,100)    
#' 
#' @export         
kinematic.viscosity <- function(Tair,pressure,constants=bigleaf.constants()){
  
  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * 1000
  
  v  <- 1.327e-05*(constants$pressure0/pressure)*(Tair/constants$Tair0)^1.81
  return(v)
}
