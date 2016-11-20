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


## Lookup/Test:
# - how to create examples
# - how to define "families"
# - how to write formulas/math. notations
# - how to write bold/italics
# - PET and typical values / how to do it right?
# - derive Ci?

# - how to deal with constants?


# ## global constants
#
constants <- list(

   cp         = 1004.834,        # specific heat of air for constant pressure [J K-1 kg-1] (Foken 2008 Eq. 2.54)
   Rgas       = 8.314510,        # universal gas constant [J mol-1 K-1] (Foken p. 245)
   Md         = 28.96,           # molar mass of dry air [g mol-1] 
   Mw         = 18.016,          # molecular mass for H2O (kg mol-1)
   eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (-) or Mw/Md
   Rv         = 461.5,           # gas constant for water vapor (J kg-1 K-1) (Stull_1988 p.641)
   Rd         = 287.0586,        # gas constant for dry air (J kg-1 K-1) (Foken p. 245)
   Kelvin     = 273.15,          # conversion Kelvin to Celsius
   gv         = 9.81,            # gravitational acceleration (m s-2)
   pressure0  = 101325,          # reference atmospheric pressure (at sea level) (Pa)
   Tair0      = 273.15,          # reference air temperature (K)
   k          = 0.41             # von Karman constant (-)

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
#' @seealso \code{\link{Gb.Su}}, \code{\link{Gb.McNaughton}}
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
#' @param wind      wind speed (m s-1)
#' @param presssure atmospheric pressure (Pa)
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
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.McNaughton}}
#' 
#' @export
Gb.Su <- function(ustar,wind,pressure,Tair,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,constants){
  Tair <- Tair + 273.15
  
  if (is.null(fc)) {
    fc  <- (1-exp(-LAI/2)) 
  } 
  
  v   <- 1.327*10^-05*(101325/pressure)*(Tair/273.15)^1.81   # kinematic viscosity (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
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
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Su}}

Gb.McNaughton <- function(ustar,leafwidth,LAI,constants){
  Rb <- 130*(sqrt(leafwidth*ustar))/LAI - 1.7
  
  #kB <- k*(120/LAI*sqrt(leafwidth*ustar) - 2.5)
  
  kB  <- Rb*constants$k*ustar
  Gb  <- 1/Rb
  
  
  return(data.frame(Rb,Gb,kB))
}



#' Air density
#' 
#' @description Air density of moist air from air temperature and pressure
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (Pa)
#' @param Mair      molar mass of air (kg mol-1)
#' @param Rgas      ideal gas constant ()  
#' 
#' Foken: p/(Rl * Tv)

air.density <- function(Tair,pressure,constants){
  Tair <- Tair + constants$Kelvin
  
  rho <- constants$Mair * pressure / (constants$Rgas * Tair) / 1000
  
  return(rho)
}


#' Barometric equation
#' 
#' @description An estimate of mean pressure at some elevation as predicted by the
#'              Barometric equation
#'              
#' # airDensity <- function(temperature , pressure) {
#   # function to calculate air density
#   # input: temperature [K]
#   #        pressure [Pa] 
#   # output: air density [kg m-3]  
#   rho <- molarMassAir * pressure / ( Rgas * temperature ) / 1000
#   return(rho)
# } 


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
Monin.Obukhov.length <- function(Tair,pressure,ustar,H,constants){
  
  rho  <- air.density(Tair,pressure,constants)
  Tair <- Tair + constants$Kelvin
  MOL  <- (-rho*constants$cp*ustar^3*Tair)/(constants$k*constants$g*H)
  
  return(MOL)
}


### from Foken_2008; p.65 (suggested by Businger_1971), in the form of Högstrom 1988
stab.correction.h.Businger <- function(zeta){
  psi_h <- numeric()
  
  y <- 0.95 * ( 1 - 11.6 * zeta)^0.5
    
  # stable
  psi_h[zeta >= 0] <- -7.8 * zeta[zeta >= 0]
  # unstable
  psi_h[zeta < 0] <- 2 * log(( 1 + y[zeta < 0] ) / 2)  
  
  return(psi_h)
    
}


### Dyer_1970; Dyer_1974 (similar:Paulson 1970)
stab.correction.h.Dyer <- function(zeta){

  psi_h <- numeric()
   
  x     <- (1 - 16 * zeta)^0.5

  # stable
  psi_h[zeta >= 0] <- -5* zeta[zeta >= 0]  # from Webb_1970, used in Magnani_1998
  # unstable
  psi_h[zeta < 0] <- 2 * log((1 + x[zeta < 0])/2)

  return(psi_h)
}





#' Aerodynamic conductance
#' 
#' @description Aerodynamic conductance, including options for Boundary layer conductance
#'              and stability correction functions.
#' 
#' @param Tair              air temperature (degC)
#' @param pressure          air pressure (Pa)
#' @param wind              wind speed (m s-1)
#' @param ustar             friction velocity (m s-1)
#' @param H                 sensible heat flux (W m-2)
#' @param zr                measurement (=reference) height (m)
#' @param zh                canopy height (m)
#' @param disp              zero-plane displacement height (m)
#' @param z0m               roughness length for momentum (m)
#' @param Dl                characteristic leaf dimension (m)
#' @param N                 number of leaf sides participating in heat exchange (1 or 2)
#' @param LAI               one-sided leaf area index (m2 m-2)
#' @param fc                fractional vegetation cover (-)
#' @param Cd                foliage drag coefficient (-)
#' @param hs                roughness length of bare soil (m)
#' @param stab_correction   stability correction function
#' @param Rb_model          boundary layer resistance formulation
#' 
#' @details
#' 
#' Aerodynamic conductance (Ga) is calculated as:
#' 
#'  Ga = 1/Ra_m + 1/Rb
#'  
#'  where Ra_m is the aerodynamic resistance for momentum and Rb the canopy boundary
#'  layer resistance.
#'  
#'  The argument stab_correction determines the stability correction function used 
#'  to account for the effect of atmospheric stability on Ra_m (Ra_m is lower for unstable
#'  and higher for stable stratification). Stratification is based on a stability parameter zeta (z-d/L),
#'  where z = reference height, d the zero-plane displacement height, and L the Monin-Obukhov length.
#'  Stability correction functions calculate the coreection based on the formulations by Businger 1971 
#'  Dyer 1970, based on zeta.
#'  If "none", atmospheric stability is neglected, and Monin-Obukhov length, and the 
#'  stability parameter zeta are not calculated.
#' 
#'  Rb_model determines the boundary layer resistance ("excess resistance") model used. Thom_1972 is an 
#'  empirical formulation based on the friction velocity (ustar):
#'  
#'    formula
#'    
#'  The model by Mcnaughton and van den Hurk 1995 (McNaughton_1995), calculates Rb
#'  based on leaf width, LAI and ustar (Note that the original formulation leaf width
#'  instead of the characteristic leaf dimension (Dl) is taken):
#'   
#'     formula
#'     
#'  The Su_2001 option calculates Rb based on the physically-based Rb model by Su et tal. 2001,
#'  a simplification of the model developed by Massmann 1999:
#'  
#'     formula
#'  
#'  The models calculate the ... parameter kB-1, which related to Rb:
#'  
#'     formula
#'  
#'  Both kB and Rb are given in the output.
#'  
#' @return a dataframe with the following columns:  
#'  
#' @export
#'    


# TO DO:
# - include stability correction functions in Ga()
# - make some constants internal (e.g. Kelvin, Rgas), and the rest as default argument
#   to the function (+ find out what "." means in front of a function!)
# - include option of constant kB
# - include option to change d and/or z0m
Ga <- function(Tair,presssure,wind,ustar,constants,H,zr,zh,disp,z0m,Dl,N,LAI,fc=NULL,Cd=0.2,hs=0.01,stab_correction=c("Dyer_1970","Businger_1971","none"),
               Rb_model=c("Thom_1972","McNaughton_1995","Su_2001")) 
  {
  
  Rb_model  <- match.arg(Rb_model)
  stability <- match.arg(stability)
  
  k <- constants$k
  
  if (Rb_model == "Thom_1972"){
    
    Gb <- Gb.Thom(ustar,constants)
    
  } else if (Rb_model == "McNaughton_1995"){
    
    Gb <- Gb.McNaughton(ustar,Dl,LAI,constants)
    
  } else if (Rb_model == "Su_2001"){
    
    Gb <- Gb.Su(ustar,wind,pressure,Tair,Dl,N,fc=fc,LAI,Cd=Cd,hs=hs,constants)
  
  } 
  
  Ra_m <- wind/ustar^2
  
  
  if (stability %in% c("Businger_1971","Dyer_1970")){
    
    disp       <- 0.7 * zh
    z0m_start  <- 0.1 * zh
    
    diff <- numeric(); ustar1 <- ustar; z0m <- numeric()
    for (i in 1:95){
      ustar1[ustar < quantile(ustar,i/100,na.rm=T)] <- NA
      z0m[i]          <- summary(nls(Ra_m ~ (log((zr - disp)/z0m)/(k*ustar1)),
                                start=c(z0m=z0m_start),na.action=na.omit,
                                data=environment()))$par[1]
      
      diff[i] <- sum((Ra_m - (log((zr - disp)/z0m[i])/(k*ustar1)))^2,na.rm=T)
    }
    
#     optim(tustar,find_z0m,)
#     
#     find_z0m(tustar,Ra_m=Ra_m,zr=20,disp=15,ustar=ustar,wind=wind)
#     
#     find_z0m <- function(tustar,Ra_m,zr,disp,k=0.41,ustar,wind){
# 
#       ustar1[ustar < quantile(ustar,tustar,na.rm=T)] <- NA
#       z0m <- summary(nls(Ra_m ~ (log((zr - disp)/z0m)/(k*ustar1)),
#                      start=c(z0m=z0m_start),na.action=na.omit,
#                      data=environment()))$par[1]
#       
#       diff <- sum((Ra_m - (log((zr - disp)/z0m[i])/(k*ustar1)))^2,na.rm=T)
#       
#       return(diff)
#     }
    
    

    print(z0m)
    
    MOL   <- Monin.Obukhov.length(Tair,pressure,ustar,H,constants)
    zeta  <- (zr-disp)/MOL
    
    if (stab_correction == "Businger_1971"){
    
      psi_h <- stab.correction.h.Businger(zeta)
    
    } else if (stab_correction == "Dyer_1970"){
    
      psi_h <- stab.correction.h.Dyer(zeta)
      
    }  
    
    Ra_m  <- pmax((log((zr - disp)/z0m) - psi_h),1e-10)/(k*ustar)
    
  } else {
    
    MOL = zeta = psi_h <- rep(NA,length=length(Ra_m))
  }
  
  
  Ga_m <- 1/Ra_m
  Ra_h <- Ra_m + Gb[,"Rb"]
  Ga_h <- 1/Ra_h
   

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Rb=Gb[,"Rb"],Gb=Gb[,"Gb"],kB=Gb[,"kB"],MOL,zeta,psi_h))

  
}
  

## test data
wind <- rnorm(100,3,0.8)
ustar <- rnorm(100,0.5,0.2)
Tair <- 25
pressure <- 100000

z <- 30
d <- 15
z0m_start <- 2
k <- 0.41

ra1 <- wind/ustar^2
ra2 <- log((z -d)/z0m)/(k*ustar)

ustar1 <- ustar
ustar1[ustar1 < 0.1] <- NA

nls(ra1[] ~ log((z -d)/z0m)/(k*ustar1),start=c(z0m=z0m_start),na.action=na.exclude)

#' Latent heat of vaporization
#' 
#' @description Latent heat of vaporization as a function of air temperature.
#' 
#' @param Tair  Air temperature (deg C)
#' 
#' @details The following formula is used:
#' 
#' \deqn{lambda = (2.501 - 0.00237Tair)10^6}
#' 
#' @return latent heat of vaporization (J kg-1) 
#' 
#' @references Stull, B., 1988: An Introduction to Boundary Layer Meteorology (p.641)
#'             Kluwer Academic Publishers, Dordrecht, Netherlands
#'             
#'             Foken_2008
#' 
#' @export
LE.vaporization <- function(Tair) {
  
    k1 <- 2.501
    k2 <- 0.00237
    
    lambda <- ( k1 - k2 * Tair ) * 1e+06
    
    return(lambda)
  }



#' Saturation vapor pressure (Esat) and slope of the Esat curve
#' 
#' @description Calculates Saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' @param Tair     air temperature (deg C)
#' @param formula  formula to be used. Either "Alduchov_1996" or "Sonntag_1990"
#' 
#' @details Esat (Pa) is calculated based on the Magnus equation:
#' 
#'  \deqn{Esat = a * exp((b * Tair) / (c + Tair))}
#'  
#'  where the coefficients a, b, c differ with the formula.
#'  The slope of the Esat curve is calculated as the first derivative of the function.
#' 
#' @return A dataframe with the following columns:
#' \itemize{
#'   \item{Esat} { - Saturation vapor pressure (Pa)}
#'   \item{delta}{ - Slope of the saturation vapor pressure curve (Pa K-1)}
#' }
#'    
#'    
#' 
#' @references Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure. Journal of Applied Meteorology, 1996, 35, 601-609
#' 
#'             Sonntag_1990 Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90, and psychrometric formulae 
#'             (--> as used by WMO 2008!)
#' @export
Esat <- function(Tair,formula=c("Alduchov_1996","Sonntag_1990")){
  
  formula <- match.arg(formula)
  
  if (formula == "Alduchov_1996"){
    a <- 610.94
    b <- 17.625
    c <- 243.04
  } else if (formula == "Sonntag_1990"){
    a <- 611.2
    b <- 17.62
    c <- 243.12
  }

  Esat <- a * exp((b * Tair) / (c + Tair))
  
  delta <- eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  
  return(data.frame(Esat,delta))
  
}



#' Conversions between humidity measures
#' 
#' @description Conversion between VPD, specific humidity, and relative humidity
#' 
#' @param e    vapor pressure (Pa)
#' @param p    air pressure (Pa)
#' @param q    specific humidity (kg kg-1)
#' @param VPD  vapore pressure deficit (kPa)
#' 
#' @family humidity conversion
e.to.q <- function(e,pressure,constants){
  q <- constants$eps * e / (pressure - 0.378 * e) 
  
  return(q)
}

#' @family humidity conversion
q.to.e <- function(q,pressure){
   e <- 5 +5
}

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




#' Psychrometric constant
#' 
#' @description Calculates the psychrometric 'constant'
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (Pa)
#' 
#' @details The psychrometric constant is given as:
#' 
#'  formula
#'  
#' @return the psychrometric constant (Pa K-1)
#'  
#' @references test
#' 
psychrometric.constant <- function(Tair,pressure,constants){
  
  lambda <- LE.vaporization(Tair)
  gamma  <- (constants$cp * pressure) / (constants$eps * lambda)
  
  return(gamma)
}

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

#' Surface conductance for water from the Penman-Monteith equation
#' 
#' @description Calculates surface conductance for H2O from the inverted Penman-Monteith
#'              equation.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (Pa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2)
#' @param S         Sum of all storage fluxes (W m-2)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (Pa)
#' @param Ga        Aerodynamic conductance (m s-1)
#' 
#' @details Surface conductance is calculated from the inverted Penman-Monteith equation:
#' 
#'  formula
#'  
#'  Note that Gs > Gc (canopy conductance) in time periods when a significant fraction of 
#'  ET comes from interception or soil evaporation. 
#'  Available energy A is given by (Rn - G - S). If G and/or S are not provided, then A = Rn. 
#' 
#' 
#' @return a dataframe with two columns: 
#' 
#'  Gs_ms:  Surface conductance in ms-1
#'  Gs_mol: Surface conductance in mol m-2 s-1
#' 
#' @export
GsPM <- function(Tair,pressure,Rn,VPD,LE,Ga,constants){

   delta <- esatFromTemperature(temperature)[,"desatdT"]
   gamma <- psychrometricConstant(temperature, pressure)
   rho   <- air.density(Tair, pressure)
 
   Gs_ms  <- ( LE * Ga * gamma ) / ( delta * Rn + rho * cp * Ga * VPD - LE * ( delta + gamma ) )

   Gs_mol <- ms.to.mol(Gs_ms,Tair,pressure,constants) 

   return(data.frame(Gs_ms,Gs_mol))

}


#' Conversion between conductance units
#' 
#' @description Converts canopy/surface conductance from ms-1
#'              to mol m-2 s-1 (ms.to.mol), or vice versa (mol.to.ms) based on
#'              the ideal gas law.
#' 
#' @param Gc_ms       Canopy/Surface conductance (ms-1)
#' @param Tair        Air temperature (deg C)
#' @param pressure    Air pressure (Pa)
#' @param constants 
#' 
#' @return Canopy/Surface conductance in mol m-2 s-1 
#' 
#' 
ms.to.mol <- function(Gc_ms,Tair,pressure,constants){
  Tair   <- Tair + constants$Kelvin
  Gc_mol <- Gc_ms * pressure / (constants$Rgas * Tair)
  
  return(Gc_mol)
}


mol.to.ms <- function(Gc_mol,Tair,pressure,constants){
  Tair  <- Tair + constants$Kelvin
  Gc_ms <- Gc_mol * (constants$Rgas * Tair) / pressure
  
  return(Gc_ms)
}


#' Conversion between LE and ET
#' 
#' @description converts water loss from mass (ET) to energy (ET) units, and vice versa
#' 
#' @param LE   Latent heat flux (W m-2)
#' @param Tair Air temperature (deg C)
#' 
#' @return ET Evapotranspiration (kg m-2 s-1)
LE.to.ET <- function(LE,Tair){

  Tair   <- Tair + Kelvin
  lambda <- LE.Vaporization(Tair)
  ET     <- LE/lambda
   
  return(ET)
}

ET.to.LE <- function(ET,Tair){
  
  Tair   <- Tair + Kelvin
  lambda <- LE.Vaporization(Tair)
  LE     <- ET*lambda
  
  return(LE)
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
# LE.vaporization <- function(temperature) {
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

# 
# 

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
