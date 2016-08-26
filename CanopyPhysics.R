#### CanopyPhysics Routine ------------------------------------------------------------
##
## Written by S. Zaehle and J. Knauer
## Date: 21.1.2016
## Version:  beta release candidate
## History: none
## Known issues: PET_Penman gives wrong results
##
####-----------------------------------------------------------------------------------

library(devtools)
library(roxygen2)
setwd("./BigLeaf")


# ## global constants
#
constants <- list(

   cp           = 1004.834,        # specific heat of air for constant pressure [J K-1 kg-1] (Foken 2008 Eq. 2.54)
   Rgas         = 8.314510,        # universal gas constant [J mol-1 K-1] (Foken p. 245)
   Mair         = 28.96,           # molar mass of dry air [g mol-1] 
   Rdryair      = 287.0586,        # gas constant of dry air [J kg-1 K-1] (Foken p. 245)
   Kelvin       = 273.15,          # conversion Kelvin to Celsius
   eps          = 0.622,           # ratio of the molecular weight of water vapor to dry air (-)
   gv           = 9.81,            # gravitational acceleration (m s-2)
   pressure0    = 101325,          # reference atmospheric pressure (at sea level) (Pa)
   Tair0        = 273.15           # reference air temperature (K)
   k            = 0.41             # von Karman constant (-)

)


# ga_complex <- function(temperature,pressure,wspeed,ustar,fh,zr,zs,Dl){
#   # computes aerodynamic conductance according to MOST following Foken's receipe
#   # source: Liu et al. 2007 HESS (based on Thom 1975)
#   # input:  temperature: air temperature [K]
#   #         pressure: air pressure [Pa]
#   #         wspeed: wind speed [m s-1]
#   #         ustar: friction velocity [m s-1]
#   #         fh: sensible heat flux [W m-2]
#   #         zr: measurement height / reference height [m]
#   #         zs: canopy height [m]
#   # output: 
#   #         aerodynamic conductance for heat [m s-1]
#   # local variables:
#   gamma_m = gamma_h <- 16           # coefficient for stability function
#   beta_m = beta_h   <- 5            # coefficient for stability function
#   z0m   <- 0.1 * zs                 # roughness length for momentum
#   disp  <- 0.67 * zs                # displacement height
#   Cd    <- 0.2                      # 
#   fc    <- 1.0                      # canopy fraction
#   hs    <- 0.002                    #
#   
#   ## excess resistance following Su et al. 2001 (based on Massman 1999)
#   kB    <- kB_Su_2001(rk,Cd,ustar,wspeed,fc,hs,pressure,standard_air_pressure_msl,temperature,Kelvin,Dl)    # excess resistance parameter according to Thom 1972
#   #kB    <- 1.35*rk*(100*ustar)^0.33
#   z0h   <- z0m/exp(kB)                                                                                      # roughness length for heat
#   
#   
#   Ga <- numeric(length(temperature))
#   Ga[] <- NA
#   valid=(wspeed>0.0&!is.na(wspeed))&(ustar>0.05&!is.na(ustar))&!is.na(fh)&!is.na(temperature)&!is.na(pressure)
#   if(!any(valid)) return(Ga)
# 
#   # density of air
#   rho <- airDensity(temperature, pressure)
#   # Monin Obukhov length (see Foken 2008: Micrometeorology formula 2.68 p.43)
#   MOL <- (-rho*cp*ustar^3*temperature)/(rk*gv*fh)
# 
#   ## stability parameters
#   zeta   <- (zr-disp)/MOL
#   zeta0m <- z0m/MOL
#   zeta0h <- z0h/MOL
# 
#   ## stability correction functions
#   psi_m <- stabilityCorrection_m(zeta,zeta0m,gamma_m,beta_m,z0m,zr,disp)
#   psi_h <- stabilityCorrection_h(zeta,zeta0h,gamma_h,beta_h,z0h,zr,disp)
# 
#   ## aerodynamic resistance for heat (ra_h), momentum (ra_m) and excess resistance (ra_b)
#   ra_h <- 1/(rk^2*wspeed)*(log((zr - disp)/z0m) - psi_m)*(log((zr-disp)/z0h) - psi_h)
#   ra_m <- 1/(rk^2*wspeed)*(log((zr - disp)/z0m) - psi_m)^2
# 
#   # excess resistance" (ra_b) depends on the ratio z0m/z0h and is adjusted
#   # through the kB-1 number which depends on vegetation type
#   # a typical value is around 2 for forests
#   ra_b <- ra_h - ra_m
# 
#   return(1/ra_h)
# }

## Alternative version
# Gb.Thom <- function(ustar,constants){
#   Rb <- 6.2*ustar^-0.667
#   Gb <- 1/Rb
#   kB <- Rb*(constants$k*ustar)
#   
#   return(cbind(Rb,Gb,kB))
# }



#' Boundary layer conductance according to Thom 1972
#' 
#' An empirical formulation for the quasi-laminar boundary layer conductance
#' based on a simple ustar dependency.
#' 
#' @param ustar friction velocity (m s-1)
#' @param k von-Karman constant (-)
#' 
#' @return a matrix with the following columns:
#' 
#' @references Thom, A., 1972: Momentum, mass and heat exchange of vegetation. Quarterly Journal of the Royal Meteorological Society, 98, 124--134
#' 
#' @seealso \code{\link{Gb.Su}}
#' 
#' 
#' @export
Gb.Thom <- function(ustar,constants){
  # kB <- 1.35*k*(100*ustar)^0.333
  Rb <- 6.2*ustar^-0.667
  Gb <- 1/Rb
  kB <- Rb*constants$k*ustar
  
  return(data.frame(Rb,Gb,kB))
}


#' Roughness Reynolds Number
#' 
#' An empirical formulation for the quasi-laminar boundary layer conductance
#' based on ustar
#' 
#' @param ustar friction velocity (m s^{-1})
#' @param k von-Karman constant (-)
#' 
#' @return a matrix containing Gb, Rb, and kB
#' 
#' @export
ReynoldsNumber <- function(hs,ustar,pressure,Tair,constants){
  # hs = roughness height of the soil [m]
  # v = kinematic viscosity of the air
  # Tair in K!!
  v  <- 1.327e-05*(constants$pressure0/pressure)*(Tair/constants$Tair0)^1.81
  Re <- hs*ustar/v
  
  return(Re)
}




#' Boundary layer conductance according to Su et al. 2001
#' 
#' A physically based formulation for the quasi-laminar boundary layer conductance. 
#'
#' @param ustar     friction velocity (m s-1)
#' @param wind         wind speed (m s-1)
#' @param Patm      atmospheric pressure (Pa)
#' @param Tair      air temperature (degC)
#' @param Dl        leaf dimension (m)
#' @param N         number of leaf sides participating in heat exchange (1 or 2)
#' @param fc        fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI       one-sided leaf area index (-)
#' @param k         von-Karman constant
#' @param Cd        foliage drag coefficient (-)
#' @param hs        roughness height of the soil (m)
#' 
#' @return a data frame with the following columns:
#' 
#' @details If fc (fractional vegetation cover) is missing, it is estimated from LAI:
#' 
#' \eqn(fc = 1 - exp(-LAI/2))
#' 
#' 
#' @references 
#' Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of two models for estimation of the roughness height for heat transfer between the land surface and the atmosphere. Journal of Applied Meteorology, 40, 1933--1951.
#' 
#' Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using 'localized near-field' Lagrangian theory. Journal of Hydrology. 223, 27--43.
#' 
            

Gb.Su <- function(ustar,wind,Patm,Tair,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,constants){
  Tair <- Tair + 273.15
  
  if (is.null(fc)) {
    fc  <- (1-exp(-LAI/2)) 
  } 
  
  v   <- 1.327*10^-05*(101325/Patm)*(Tair/273.15)^1.81   # kinematic viscosity (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
  Re  <- hs*ustar/v                                      # Reynolds number
  kBs <- 2.46*(Re)^0.25 - log(7.4)                       # Su_2001 Eq. 13       
  Reh <- Dl*wind/v                                       # Reynolds number
  Ct  <- 1*0.71^-0.6667*Reh^-0.5*N                       # heat transfer coefficient of the leaf (Massman_1999 p.31)
  
  kB  <- (constants$k*Cd)/(4*Ct*ustar/wind)*fc^2 + kBs*(1 - fc)^2
  Rb  <- kB/(constants$k*ustar) 
  Gb  <- 1/Rb
  
  
  return(data.frame(Rb,Gb,kB))
}


#' Boundary layer conductance according to McNaughton and van den Hurk 1995 --> check units in original paper!
#' 
#' @param ustar     friction velocity (m s-1)
#' @param leafwidth leaf width (m)
#' @param LAI       one-sided leaf area index
#' @param k         von-Karman constant (-)
#' 
#' 

Gb.McNaughton <- function(ustar,leafwidth,LAI,constants){
  Rb <- 130*(sqrt(leafwidth*ustar))/LAI - 1.7
  
  #kB <- k*(120/LAI*sqrt(leafwidth*ustar) - 2.5)
  
  kB  <- Rb*constants$k*ustar
  Gb  <- 1/Rb
  
  
  return(data.frame(Rb,Gb,kB))
}



#' Air density
#' 
#' Air density from air temperature and pressure
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (Pa)
#' @param Mair      molar mass of air (kg mol-1)
#' @param Rgas      ideal gas constant ()  
#' 
#' 
#' # airDensity <- function(temperature , pressure) {
#   # function to calculate air density
#   # input: temperature [K]
#   #        pressure [Pa] 
#   # output: air density [kg m-3]  
#   rho <- molarMassAir * pressure / ( Rgas * temperature ) / 1000
#   return(rho)
# } 
AirDensity <- function(Tair,pressure,constants){
  Tair <- Tair + constants$Kelvin
  
  rho <- constants$Mair * pressure / (constants$Rgas * Tair) / 1000
  
  return(rho)
}


#' Monin-Obukhov length
#' 
#' title
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (Pa)
#' @param ustar     friction velocity (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param constants constants required
#'                  cp
#'                  k
#'                  g 
#' 
#' MOL <- (-rho*cp*ustar^3*temperature)/(rk*gv*fh)
MoninObukhovlength <- function(Tair,pressure,ustar,H,constants){
  
  rho  <- AirDensity(Tair,pressure,constants)
  Tair <- Tair + constants$Kelvin
  MOL  <- (-rho*constants$cp*ustar^3*Tair)/(constants$k*constants$g*H)
  
  return(MOL)
}


#' Aerodynamic conductance
#' 
#' Title
#' 
#' @param Tair        air temperature (degC)
#' @param pressure    air pressure (Pa)
#' @param wind        wind speed (m s-1)
#' @param ustar       friction velocity (m s-1)
#' @param H           sensible heat flux (W m-2)
#' @param zr          measurement (=reference) height (m)
#' @param zh          canopy height (m)
#' @param disp        zero-plane displacement height (m)
#' @param z0m         roughness length for momentum (m)
#' @param Dl          characteristic leaf dimension (m)
#' @param N           number of leaf sides participating in heat exchange (1 or 2)
#' @param fc          fractional vegetation cover (-)
#' @param LAI         one-sided leaf area index (m2 m-2)
#' @param k           von-Karman constant (-)
#' @param Cd          foliage drag coefficient (-)
#' @param hs          roughness length of bare soil (m)
#' @param stability   stability correction function
#' @param Rbmodel     boundary layer resistance formulation
#' 
#' @details 1. McNaughton: Dl instead of leafwidth is taken

Ga <- function(Tair,presssure,wind,ustar,H,zr,zh,disp,z0m,Dl,N,fc=NULL,LAI,k=0.41,Cd=0.2,hs=0.01,stability=c("Businger_1971","Dyer_1970","none"),
               Rbmodel=c("Thom_1972","McNaughton_1995","Su_2001"))
  
  {
  if (Rbmodel == "Thom_1972"){
    
    Rb <- 6.2*ustar^-0.667
    Gb <- 1/Rb
    kB <- Rb*k*ustar
    
  } else if (Rbmodel == "McNaughton_1995"){
    
    Rb <- 130*(sqrt(Dl*ustar))/LAI - 1.7
    Gb <- 1/Rb
    kB <- Rb*k*ustar
    
  } else if (Rbmodel == "Su_2001"){
    
    Tair <- Tair + 273.15
    
    if (is.null(fc)) {
      fc  <- (1-exp(-LAI/2)) 
    } 
    
    v   <- 1.327*10^-05*(101325/Patm)*(Tair/273.15)^1.81   # kinematic viscosity (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
    Re  <- hs*ustar/v                                      # Reynolds number
    kBs <- 2.46*(Re)^0.25 - log(7.4)                       # Su_2001 Eq. 13       
    Reh <- Dl*wind/v                                       # Reynolds number
    Ct  <- 1*0.71^-0.6667*Reh^-0.5*N                       # heat transfer coefficient of the leaf (Massman_1999 p.31)
    
    kB  <- (k*Cd)/(4*Ct*ustar/wind)*fc^2 + kBs*(1 - fc)^2
    Rb  <- kB/(k*ustar) 
    Gb  <- 1/Rb
    
  }
  
  
  z0h <- z0m/exp(kB)   ## add: if z0m not provided, it is estimated from u/ustar^2
  
  if (stability == "none"){
    
    Ra_m <- wind/ustar^2
    
  } else {
    
    
    
  } 
    
    if (stability == "Businger_1971"){
    
    psi_h <- stabilityCorrection_h3(zeta)
    Ra_m  <- pmax((log((zr - disp)/z0m) - psi_h),0.00001)/(k*ustar)
    
  } else if (stability == "Dyer_1970"){
    
    psi_h <- stabilityCorrection_h2(zeta)
    Ra_m  <- pmax((log((zr - disp)/z0m) - psi_h),0.00001)/(k*ustar)
  }  
    
    Ga_m <- 1/Ra_m
    Ra_h <- Ra_m + Rb
    Ga_h <- 1/Ra_h
    
  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Gb,Rb,kB))
    
  }
  
  
  
  
  
  
  }


# ga_complex <- function(temperature,pressure,wspeed,ustar,fh,zr,zs,Dl,disp,z0m,LAI,kB_version,stab_version){
#   # computes aerodynamic conductance according to MOST following Foken's receipe
#   # source: Liu et al. 2007 HESS (based on Thom 1975)
#   # input:  temperature: air temperature [K]
#   #         pressure: air pressure [Pa]
#   #         wspeed: wind speed [m s-1]
#   #         ustar: friction velocity [m s-1]
#   #         fh: sensible heat flux [W m-2]
#   #         zr: measurement height / reference height [m]
#   #         zs: canopy height [m]
#   #         Dl: leaf dimension [m]
#   #         disp: displacement height [m]
#   #         z0m: roughness length for momentum [m]
#   #         kB_version: model used to calculated kB (excess resistance parameter kB): either "Su_2001" or "Thom_1972"
#   #         stab_version: stability function used: "Dyer_1970" or "Businger_1971" or "none"
#   # output: 
#   #         aerodynamic conductance for heat [m s-1]
#   # local variables:
#   Cd    <- 0.2                      # 
#   hs    <- 0.002                    # roughness length of bare soil 
#   
#   
#   ## excess resistance
#   if (kB_version == "Su_2001"){
#     kB  <- kB_Su_2001(rk,Cd,ustar,wspeed,LAI,hs,pressure,standard_air_pressure_msl,temperature,Kelvin,Dl) # Su et al. 2001 (based on Massman 1999)
#   } else if (kB_version == "Thom_1972"){
#     kB <- 1.35*rk*(100*ustar)^0.33
#   }  
#   z0h   <- z0m/exp(kB) # roughness length for heat
# 
#   # density of air
#   rho <- airDensity(temperature, pressure)
#   # Monin Obukhov length (see Foken 2008: Micrometeorology formula 2.68 p.43)
#   MOL <- (-rho*cp*ustar^3*temperature)/(rk*gv*fh)
# 
#   ## stability parameters
#   zeta   <- (zr-disp)/MOL
#   zeta0m <- z0m/MOL
#   zeta0h <- z0h/MOL
# 
#   ## stability correction functions
#   #psi_m <- stabilityCorrection_m(zeta,zeta0m,gamma_m,beta_m,z0m,zr,disp)
#   #psi_h <- stabilityCorrection_h(zeta,zeta0h,gamma_h,beta_h,z0h,zr,disp)
#   
#   if (stab_version == "Dyer_1970"){
#     psi_h <- stabilityCorrection_h2(zeta)
#   } else if (stab_version == "Businger_1971"){
#     psi_h <- stabilityCorrection_h3(zeta)
#   } else if (stab_version == "none"){
#     psi_h <- rep(0,length(zeta))
#   }
#   
#   
#   ## aerodynamic resistance for heat (ra_h), momentum (ra_m) and excess resistance (ra_b)
#   ra_h  <- pmax((log((zr - disp)/z0m) - psi_h),0.00001)/(rk*ustar) + kB/(rk*ustar)      # from Verma_1989
#   ra_m  <- pmax((log((zr - disp)/z0m) - psi_h),0.00001)/(rk*ustar)
#   #ra_m <- 1/(rk^2*wspeed)*(log((zr - disp)/z0m) - psi_m)^2
# 
#   ra_b <- ra_h - ra_m
# 
#   return(cbind(1/ra_h,1/ra_m,1/ra_b,kB,z0h,psi_h))
# }
# 
# 
# 
# ga_simple <- function(wspeed,ustar){
#   # as previously calculated in the database, neglecting stability corrections
#   # second term represents excess resistance (from Monteith & Unsworth 2010 p. 341 Formel 17.8)
#   ra_m <- wspeed/ustar^2
#   ra_b <- 6.2*ustar^-0.667
#   ra_h <- ra_m + ra_b
#   ga <- 1/ra_h
#   return(ga)
# }
# 
# 
# ga_simple_ustar <- function(wspeed,ustar,zr,z0m,disp){
#   #kB  <- 1.35*rk*(100*ustar)^0.33
#   #z0h <- z0m/exp(kB)      # roughness length for heat
# 
#   ra <- log((zr - disp)/z0m)/(rk*ustar)
#   ga <- 1/ra
#   return(ga)
# }
# 
# 
# 
# gcPenmanMonteithMethod <- function(temperature,pressure,rnet,vpd,fle,ga){
#   # function to calculate surface conductance by inverting the Penman Moneith equation
#   # source: Monteith & Unsworth, 2008
#   # input:  temperature: temperature [K]
#   #         pressure: pressure [Pa]
#   #         rnet: net radiation [W m-2]
#   #         vpd: vapour pressure deficit [Pa]
#   #         fle: latent heat flux [W m-2]
#   #         ga: aerodynamic conductance [m s-1]
#   # output  gc: surface conductance [m s-1]
#   #         gc_mol: surface conductance [mol s-1]  
#   valid = (fle>0&!is.na(fle)&rnet>0&!is.na(rnet)&fle<rnet)
# 
#   delta <- esatFromTemperature(temperature)[,"desatdT"]
#   gamma <- psychrometricConstant(temperature, pressure)
#   rho <- airDensity(temperature, pressure)
# 
#   gc <- ( fle * ga * gamma ) / ( delta * rnet + rho * cp * ga * vpd - fle * ( delta + gamma ) )
#   gc[!valid] <- NA
#   gc[gc>ga|gc<=0] <- NA
#   gc_mol <- convertToMol(gc,temperature,pressure) 
# 
#   return(cbind(gc,gc_mol))
# 
# }
# 
# 
# 
# petPenmanMethod <- function(temperature, pressure, rnet, vpd, ga) {
#   # function to calculate PET given Ga
#   # source: Monteith & Unsworth, 2008, eq 13.28
#   # input:  temperature [K]
#   #         pressure [Pa]
#   #         rnet: net radiation [W m-2]
#   #         vpd: vapour pressure deficit [Pa]
#   #         ga: aerodynamic conductance for water [m s-1]  
#   # output  pet: instantaneous potential evapotranpiration [mm s-1] 
#   delta <- esatFromTemperature(temperature)[,"desatdT"] # [Pa K-1]
#   gamma <- psychrometricConstant(temperature, pressure) # [Pa K-1]
#   lambda <- latentHeatOfVapourisation(temperature)      # [J kg-1]
#   rho <- airDensity(temperature, pressure)              # [kg m-3]
#   pet <- ( delta * rnet + rho * cp * vpd * ga ) / (( delta + 0.93 * gamma ) *  lambda )
#   return(pet)
# }
# 
# petPriestleyTaylorMethod <- function(temperature, pressure, rnet) {
#   # function to calculate PET given Ga
#   # source: Monteith & Unsworth, 2008, eq 13.38
#   # input:  temperature [K]
#   #         pressure [Pa]
#   #         rnet: net radiation [W m-2]
#   # output  pet: instantaneous potential evapotranpiration [mm s-1] 
#   # local constant:
#   k_pt = 1.26
#   # derived parameters:
#   delta <- esatFromTemperature(temperature)[,"desatdT"] # [Pa K-1]
#   gamma <- psychrometricConstant(temperature, pressure) # [Pa K-1]
#   lambda <- latentHeatOfVapourisation(temperature)      # [J kg-1]
# 
#   # Priestly Taylor equation
#   pet <- k_pt * delta * rnet / (delta + gamma) / lambda 
#   return(pet)
# }
# 
# 
# calcSurfaceConditions <- function(temperature, pressure, rnet, fh, fle, vpd, ga) {
# 
#   rho <- airDensity(temperature, pressure)
#   # inferred surface temperature from energy budget
#   surface_temperature <- temperature + fh / ( rho * cp * ga)
#   #surface_temperature[abs(surface_temperature-temperature)>10] <- NA
# 
#   # penman potential evaporation
#   pet_penman <- petPenmanMethod(temperature, pressure, rnet, vpd, ga)
#   pet_pt <- petPriestleyTaylorMethod(temperature, pressure, rnet)
# 
#   # vapour pressure deficit at canopy level
#   # solving the equation E = ga*rho * (qair - qsurf) in mass notation or pressure
#   gamma <- psychrometricConstant(temperature,pressure)
#   
#   esat     <- esatFromTemperature(temperature)[,'esat']
#   esatsurf <- esatFromTemperature(surface_temperature)[,'esat']
#   
#   eair  <- esat - vpd
#   esurf <- eair + (fle * gamma)/(ga * rho * cp)
#   
#   surface_vpd <- esatsurf - esurf
#   surface_vpd[surface_vpd<0|surface_vpd>esatsurf]<-NA 
#   
#   ## calculate specific humidities
#   qsurf    <- convertVapourPressureToSpecHum(esurf,pressure)
#   qsatsurf <- convertVapourPressureToSpecHum(esatsurf,pressure)
# 
#   return(cbind(surface_temperature,esurf,esatsurf,qsurf,qsatsurf,surface_vpd,pet_penman,pet_pt))
# }
# 
# 
# ### helper functions ------------------------------------------------###
# 
# ## Foken Micrometeorology formula 2.69
# # Temp: temperature       [K]
# #    q: specific humidity [kg kg-1]
# Temp_virtual <- function(Temp,q){
#   Tv <- Temp*(1+0.61*q)   # mixing ratio is approximated by specific humidity
#   return(Tv)
# }
# 
# airDensity <- function(temperature , pressure) {
#   # function to calculate air density
#   # input: temperature [K]
#   #        pressure [Pa] 
#   # output: air density [kg m-3]  
#   rho <- molarMassAir * pressure / ( Rgas * temperature ) / 1000
#   return(rho)
# } 
# 
# convertToMol <- function(conductance,temperature,pressure) {
#   # function to convert conductance from m s-1 to mol m-2 s-1
#   # input:  conductance [m s-1]
#   #         temperature [K]
#   #         pressure [Pa]  
#   conductance_mol <- conductance * pressure / ( Rgas * temperature ) 
#   return(conductance_mol)
# }
# 
# convertSpecHumToVapourPressure <- function(spec_hum,pressure) {
#   # function to convert specific humidity into vapour pressure
#   # source: Monteith & Unsworth 2008
#   # input:  spec_hum: specific humidity [g g-1]
#   #         pressure: atmospheric pressure [Pa]
#   # output: vapour_pressure: vapour pressure [Pa] 
#   #vapour_pressure = spec_hum * pressure / ( spec_hum + eps * ( 1 + spec_hum ))
#   vapour_pressure <- spec_hum / eps * pressure
#   return(vapour_pressure)
# }
# 
# convertVapourPressureToSpecHum <- function(vapour_pressure,pressure) {
#   # function to convert specific humidity into vapour pressure
#   # source: Monteith & Unsworth 2008
#   # input:  vapour_pressure: vapour pressure [Pa]
#   #         pressure: atmospheric pressure [Pa]
#   # output: spec_hum: specific humidity [g g-1]
#   spec_hum <- eps * vapour_pressure / pressure 
#   return(spec_hum)
# }
# 
# esatFromTemperature <- function(temperature) {
#   # function to calculate saturating vapour pressure for a given temperature
#   # source: Monteith & Unworth 1995. Eq. 2.27
#   # input: temperature [K]
#   # output: saturating vapour pressure [Pa]
#   #         slope of the saturating vapour pressure [Pa K-1]
#   # local constants
#   eps <- 611   
#   A <- 17.27  
#   B <- 237.15 
#   d0 <- 4089 
#   # saturating vapour pressure 
#   esat <- eps * exp (A * (temperature - Kelvin) / (temperature - Kelvin + B ))
#   # slope of the saturation vapour pressure curve w.r.t. temperature
#   desatdT <- d0 * esat / ((  temperature - Kelvin + B ) ^ 2 ) 
#   return(cbind(esat,desatdT))
# }
# 
# latentHeatOfVapourisation <- function(temperature) {
#   # function to calculate the latent heat of vapourisation
#   # source Stull, 1988, p 641
#   # Note: this version give intermediate results to the tables of Monteith/Unsworth 2008,
#   #       and Bonan 2008
#   # input:  temperature [K]
#   # output: latent heat of vapourisation [J kg-1]
#   # local constants
#   k1 <- 2.501
#   k2 <- 0.00237
#   lambda <- ( k1 - k2 * ( temperature - Kelvin )) * 1e+06
#   return(lambda)
# }
# 
# decouplingFactor <- function(temperature,pressure,ga,gc) {
#   # function to calculate omega, the decoupling factor of vegetation
#   # source: Jarvis & McNaughton (1986): Stomatal Control of Transpiration:
#   #         scaling up from leaf to region (Advances in Ecological Research,15) 
#   # input:  temperature: temperature [K]
#   #         pressure: pressure [Pa]
#   #         ga: aerodynamic conductance [m s-1]
#   #         gc: canopy / surface conductance [m s-1] 
#   # output: omega [dimension less]
#   delta <- esatFromTemperature(temperature)[,"desatdT"]
#   gamma <- psychrometricConstant(temperature,pressure)
#   epsilon <- delta/gamma   
#   omega <- (epsilon + 1)/(epsilon + 1 + ga/gc)
#   return(omega)
# }
# 
# 
# ## JK: using the hypsometric equation would be better here, but only
# ##     if specific humidity is available from the EC stations
# pressureFromAltitude <- function(altitude) {
#   # function to calculate air pressure based on altitude
#   # source: Stull (2000): Meteorology for Scientists and Engineers, p.12,13
#   # input: altitude [m]
#   # output: atmospheric pressure [Pa]
#   # local constants
#   hlp1 <- 288.15
#   hlp2 <- 6.5
#   hlp3 <- 1000
#   hlp4 <- -5.255877
#   # pressure from altitude
#   dum <- hlp1 - (hlp2 * altitude / hlp3 )
#   air_pressure <- standard_air_pressure_msl * (hlp1 / dum) ^ hlp4
#   return(air_pressure)
# }
# 
# 
# barometric_formula <- function(Altitude,Temp,VPD){
#   # calculate surface pressure with the barometric formula
#   # first approximate pressure with air temperature,
#   # then use virtual temperature to calculate final pressure
#   # input: altitude [m]
#   #        air temperature [K]
#   #        vapor pressure deficit [Pa]
#   # output: atmospheric pressure at surface [Pa]
#   p0 <- standard_air_pressure_msl
#   
#   p  <- p0 * exp((-gv*Altitude)/(R_dryair * Temp))
#   
#   Esat <- as.numeric(esatFromTemperature(Temp)[,1])
#   eact <- Esat - VPD
#   q    <- convertVapourPressureToSpecHum(vapour_pressure=eact,pressure=p)
#   
#   Tv <- Temp_virtual(Temp,q)
#   
#   p  <- p0 * exp((-gv*Altitude)/(R_dryair * Tv))
#   return(p)
# }
# 
# 
# psychrometricConstant <- function(temperature, pressure) {
#   # function to calculate the psychrometric constant
#   # source: Bonan et al. 2008
#   # input:  temperature [K]
#   #         pressure [Pa]
#   # output: psychrometric constant gamma [Pa K-1]
#   lambda <- latentHeatOfVapourisation(temperature)
#   gamma <- ( cp * pressure ) / (eps * lambda)
#   return(gamma)
# }
# 
# stabilityCorrection_h <- function(zeta,zeta0h,gamma_h,beta_h,z0h,z,d){
#   # computes the stability function for momentum according to the MOST
#   # source: Liu et al. 2007 HESS, p 771 
#   psi_h <- as.numeric(rep(NA,length(zeta)))
# 
#   # helping variables
#   y  <- (1 - gamma_h*zeta)^0.5
#   y0 <- (1 - gamma_h*zeta*(z0h/(z-d)))^0.5
# 
#   for (i in seq_along(psi_h)){
#     if(!is.na(zeta[i])) { 
#       if (zeta[i] >= 0){    # stable conditions (zeta >= 0)
#          psi_h[i] <- -beta_h*(zeta[i] - zeta0h[i])
#       } else {              # unstable conditions (zeta < 0)
#          psi_h[i] <- 2*log((1 + y[i])/(1 + y0[i]))
#       }
#     }
#   }
#   return(psi_h)
# }
# 
# stabilityCorrection_h2 <- function(zeta){
#   ## Dyer_1970; Dyer_1974 (similar:Paulson 1970)
#   psi_h <- as.numeric(rep(NA,length(zeta)))
#   
#   x     <- (1 - 16*zeta)^0.5
#   
#   for (i in seq_along(psi_h)){
#     if(!is.na(zeta[i])) { 
#       if (zeta[i] >= 0){    # stable conditions (zeta >= 0) from Webb_1970, used in Magnani_1998
#         psi_h[i] <- -5*zeta[i]
#       } else {              # unstable conditions (zeta < 0) from Dyer_1970
#         psi_h[i] <- 2*log((1 + x[i])/2)
#       }
#     }
#   }
#   return(psi_h)
#   
# }
# 
# 
# stabilityCorrection_h3 <- function(zeta){
#   ## from Foken_2008; p.65 (suggested by Businger_1971)
#   psi_h <- as.numeric(rep(NA,length(zeta)))
#   
#   #x <- (1- 19.3*zeta)^0.25
#   y <- 0.95*(1-11.6*zeta)^0.5
#   
#   for(i in seq_along(psi_h)){
#     if(!is.na(zeta[i])){
#       if(zeta[i] >= 0){   # stable conditions (zeta >= 0)
#         psi_h[i] <- -7.8*zeta[i]    
#       } else {            # unstable conditions
#         psi_h[i] <- 2*log((1+y[i])/2)
#       }
#     }
#   }
#   return(psi_h)
# }
# 
# #one   <- (1-16*zeta)^0.5         # Magnani_1998, Dyer_1974
# #two   <- 0.95*(1-11.6*zeta)^0.5  # Foken_2008
# #three <- ((1 - 16*zeta)^0.25)^2    # Paulson_1970
# 
# 
# stabilityCorrection_m <- function(zeta,zeta0m,gamma_m,beta_m,z0m,z,d){
#   # computes the stability function for momentum according to the MOST
#   # source: Liu et al. 2007 HESS, p771 
#   psi_m <- as.numeric(rep(NA,length(zeta)))
# 
#   # helping variables
#   x  <- (1 - gamma_m*zeta)^0.25
#   x0 <- (1 - gamma_m*zeta*(z0m/(z-d)))^0.25
# 
#   for (i in seq_along(psi_m)){
#     if(!is.na(zeta[i])) { 
#        if (zeta[i] >= 0){    # stable conditions (zeta >= 0)
#          psi_m[i] <- -beta_m*(zeta[i] - zeta0m[i])
#        } else {              # unstable conditions (zeta < 0)
#          psi_m[i] <- 2*log((1 + x[i])/(1 + x0[i])) + log((1 + x[i]^2)/(1 + x0[i]^2)) - 2*atan(x[i]) +
#            2*atan(x0[i])
#        }
#     }
#   }
#   return(psi_m)
# }
# 
# 
# 
# 
# ### Roughness Reynolds number
# Reynolds_number <- function(hs,ustar,p,p0=p0,Tair,Tair0=Kelvin){
#   # hs = roughness height of the soil [m]
#   # v = kinematic viscosity of the air
#   # Tair in K!!
#   v  <- 1.327e-05*(p0/p)*(Tair/Tair0)^1.81
#   Re <- hs*ustar/v
#   return(Re)
# }
# 
# 
# ### kB model from Su_2001
# kB_Su_2001 <- function(rk,Cd,ustar,u,LAI,hs,p,p0,Tair,Tair0,Dl){
#   fc  <- (1-exp(-LAI/2))                                                        ## as in JSBACH
#   Re  <- Reynolds_number(hs=hs,ustar=ustar,p=p,p0=p0,Tair=Tair,Tair0=Kelvin)    ## Su 2001 according to Brutsaert_1982
#   kBs <- 2.46*(Re)^0.25 - log(7.4)                                              # Su_2001 Eq. 13
#   
#   v   <- 1.327*10^-05*((p0/1000)/(p/1000))*(Tair/Tair0)^1.81                    # kinematic viscosity (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
#   Reh <- Dl*u/v                                                                 # Reh = Reynolds number
#   Ct  <- 1*0.71^-0.6667*Reh^-0.5*2                                              # heat transfer coefficient of the leaf (Massman_1999 p.31)
#   kB  <- (rk*Cd)/(4*Ct*ustar/u)*fc^2 + kBs*(1 - fc)^2
#   return(kB)
# }
# 
# 
# # 
# # kinem_viscosity <- function(p0=101325,p=100000,Tair=26.85,Tair0=273.15){
# #   v <- 1.327*10^-05*((p0/1000)/(p/1000))*((Tair+Kelvin)/Tair0)^1.81
# #   return(v)
# # }
# # 
# # kinem_viscosity(p0=101300,p=101300)
# 
# # # 
# # u     <- rnorm(1,2,0.2)
# # ustar <- rnorm(1,0.5,0.05)
# # Dl    <- 0.06
# # Tair  <- 20
# # LAI   <- 4
# # hs    <- 0.005
# 
# # kB_Su_2001(rk=0.41,Cd=0.2,ustar,u,LAI=1,hs=0.005,100000,101325,Tair=293,Tair0=273,Dl=0.05)
# 
# 
# 
# ### compare stability functions
# # zeta <- seq(-4,1,0.01)
# # zeta0h <- 0.2*zeta
# # psi1 <- stabilityCorrection_h(zeta,zeta0h=zeta0h,gamma_h=16,beta_h=5,z0h=2,z=20,d=14)
# # psi2 <- stabilityCorrection_h2(zeta)
# # psi3 <- stabilityCorrection_h3(zeta)
# # 
# # plot(psi2 ~ zeta)
# # points(psi1 ~ zeta,col="blue")
# # points(psi3 ~ zeta,col="red")
# 
# 
# LEtoET <- function(LE,temperature){
#   # Input:
#   # latent heat flux LE [W m-2]
#   # air temperature [K]
#   
#   # Output
#   # ET in kg m-2 s-1
#   
#   lambda <- latentHeatOfVapourisation(temperature)
#   ET     <- LE/lambda
#   
#   return(ET)
# }
# 
# 
# 


###################################################################################
#### Energy balance ###############################################################
###################################################################################

#' biochemical energy (Heat)
#' 
#' radiant energy absorbed in photosynthesis or heat release by respiration based on NEE  
#' 
#' @param alpha   energy taken up/released by photosynthesis/respiration (J umol-1)
#' @param NEE     net ecosystem exchange (umol CO2 m-2 s-1)
#' 
#' @details the following sign convention is employed: NEE is negative when carbon is taken up by the ecosystem.
#'          Positive values in the result mean that energy is taken up by the ecosystem, negative ones that heat is 
#'          released.
#'          The value of alpha is taken from Nobel_1974 (Meyers_2004), but other values exist (Blanken_1997)
#' 
#' @return biochemical energy Sp (W m-2)

Photosyn_energy <- function(alpha=0.422,NEE){
  Sp <- -alpha*NEE
  
  return(Sp)
}








## General literature
## Gb: Hong_2012


# 
