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
# - how to write bold/italics
# - PET and typical values / how to do it right?
# - derive Ci? 
# - include extrapolation of Gs?
# - would it be better to convert data.frame right at the beginning and then 
#   save a few lines of code where columns are checked and converted to vectors?
#   or alternatively: combine to a smaller number of functions, so that this
#   check does not have to be performed so many times? See how Twutz did it!
#
# - include aggregation in data or enable implementation in wrapper function 
#   such as aggregate()

### Issues/check:
# MOL independent of temperature?
# more efficient way to provide data.frames? rather than specifiyng each column separately?
# check if Martin 1989, Daudet 1999 etc. really need Gb only or Ga?
# read + understand the differences of the omega formulations, starting with appendix of Jarvis & McNaughton 1986!

## TODO:
# Esat <- provide only one option? Else this option has to be inlcuded in all functions where Esat is used!!
#         not really, if it was  to be changed, do it in the function directly, not where it is a sub-function!
# check e to q conversions again!
# References for EVERY function!!!
# check for more functions in Jarvis,McNaughton and the Monteith textbook!!!



### still missing:
# Ga complete 
# the estimation of d and z0m; reference with dependency on LAI? but not necessary now!
# wind speed at the top of the canopy
# ETpot - add other formulations: see Donahue 2010!
# ET_imp, ET_equ
# Ci
# inst. energy balance closure + correction
# Uncertainty!!





## test data
wind     <- rnorm(20,3,0.8)
ustar    <- rnorm(20,0.5,0.2)
Tair     <- rnorm(20,25,2)
pressure <- rep(101.325,20)
Rn       <- rnorm(20,400,40)
H        <- 0.3*Rn
LE       <- 0.4*Rn
G        <- 0.1*Rn
Ga       <- rnorm(20,0.08,0.01)
Gs       <- rnorm(20,0.3,0.02)*(1/41)
VPD      <- rnorm(20,2,0.3)
precip   <- c(rep(0,18),0.1,0.001)
PPFD     <- seq(180,300,length.out=20)
GPP      <- rnorm(20,20,3)
NEE      <- -GPP
day      <- c(rep(1,10),rep(2,10))



testdf <- data.frame(wind,ustar,Tair,pressure,Rn,
                     H,LE,G,Ga,Gs,VPD,precip,PPFD,GPP,NEE,day)





#' checks columns in a data frame 
#' 
#' @description check if columns in a data.frame or matrix exist and if they 
#'              are of the right type (numeric).
#' 
#' @param data         a data.frame or a matrix
#' @param column_name  the column name
#' 
#' 
#' 
check.columns <- function(data,column_name){
  if (column_name %in% colnames(data)){
    var <- data[,column_name]
    if (is.numeric(var)){
      return(var)
    } else {
      stop("variable '",column_name,"' has to be numeric")
    }
  } else {
    stop("column '",column_name,"' does not exist")
  }
}


# ## global constants
bigleaf.constants <- function(){
  
  list(
    cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1) (Foken 2008 Eq. 2.54)
    Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1) (Foken p. 245)
    Md         = 0.0289645,       # molar mass of dry air [kg mol-1] (Foken 2008) 
    Mw         = 0.0180153,       # molar mass of water vapor [kg mol-1] (Foken 2008)
    eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (-) or Mw/Md
    Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull_1988 p.641)
    Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken p. 245)
    Kelvin     = 273.15,          # conversion degree Celsius to Kelvin
    g          = 9.81,            # gravitational acceleration (m s-2)
    pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
    Tair0      = 273.15,          # reference air temperature (K)
    k          = 0.41,            # von Karman constant (-)
    Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
    Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
    sigma      = 5.670367e-08     # Stefan-Boltzmann constant (W m-2 K-4)
  )

}


########################
### Filter functions ###----------------------------------------------------------------------------
########################

## TODO:
# possibility to provide time information as POSIX or similar
# warning if neither daily nor hourly
# include option: number of NEW data excluded!

#' Filter data
#'
#' @description Filters timeseries in which Gs describes physiologically meaningful
#'              quantity
#' 
#' @param data          data.frame or matrix containing all required input variables in 
#'                      half-hourly or hourly resolution. Including year, month, day information
#' @param precip        precipitation (mm)
#' @param PPFD          photosynthetically photon flux density (umol m-2 s-1)
#' @param Tair          air temperature (deg C)
#' @param ustar         friction velocity (m s-1)
#' @param VPD           vapor pressure deficit (kPa)
#' @param rH            relative humidity (-), only one of VPD and rH required
#' @param GPP           gross primary productivity (umol m-2 s-1)
#' @param tprecip       precipitation threshold used to demark a precipitation event (mm)
#' @param precip_hours  number of hours removed following a precipitation event (h)
#' @param trad          radiation threshold (umol m-2 s-1)
#' @param ttemp         temperature threshold (deg C)
#' @param tustar        friction velocity threshold (m s-1)
#' @param tGPP          GPP threshold (???)
#' @param trH           relative humidity threshold (-)
#' @param NAasprecip    if TRUE missing precipitation data are treated as precipitation events
#' 
#' 
#' @details growing season:
#'          growing season is filtered based on daily GPP (aggregated from halfhourly values).
#'          
#' @return filtered_data the same data frame as provided as argument to the function
#'                       but with filtered timeperiods set to NA in all columns except
#'                       the ones used for filtering. The data.frame contains one additional
#'                       column "valid", which indicates whether the timestep is valid 
#'                       according to the filter criteria (1) or invalid (0).
#' 
#' @references Knauer 2017
#' 
#' @export                     

filter.data <- function(data,precip="precip",PPFD="PPFD",Tair="Tair",ustar="ustar",VPD="VPD",GPP="GPP",day="day",month="month",year="year",
                        tprecip=0.01,precip_hours=24,trad=250,ttemp=5,tustar=0.2,tGPP=0.5,trH=0.95,NAasprecip=F){
  
  ## test data availability
  precip <- check.columns(data,precip)
  PPFD   <- check.columns(data,PPFD)
  Tair   <- check.columns(data,Tair)
  ustar  <- check.columns(data,ustar)
  VPD    <- check.columns(data,VPD)
  GPP    <- check.columns(data,GPP)
  day    <- check.columns(data,day)
  month  <- check.columns(data,month)
  year   <- check.columns(data,year)
  
  rH <- VPD.to.rH(VPD,Tair)
  
  ## hourly or halfhourly data?
  tstep_day <- as.integer(nrow(data) / length(unique(day)))
  thour     <- ifelse(tstep_day == 24,1,2)
  
  valid <- rep(1L,nrow(data))
  
  ## start filtering
  # 1) GPP
  GPP_daily <- aggregate(GPP,by=list(day),mean,na.rm=T)[,2]
  growseas_invalid <- filter.growseas(GPP_daily,tGPP=tGPP,ws=12,min_int=3)
  growseas_invalid <- which(day %in% growseas_invalid)
  
  # 2) precipitation
  if (NAasprecip){
    precip_events <- which(precip > tprecip | is.na(precip))
  } else {
    precip_events <- which(precip > tprecip)
  }
  precip_invalid <- unique(as.numeric(unlist(sapply(precip_events, function(x) x:(min(x+precip_hours*thour,nrow(data),na.rm=T))))))
  
  # 3) meteorological variables (PPFD, Tair, ustar, rH) 
  PPFD_invalid  <- which(PPFD <= trad)
  Tair_invalid  <- which(Tair <= ttemp)
  ustar_invalid <- which(ustar <= tustar)
  rH_invalid    <- which(rH >= trH)
  
  # 4) count number of filtered data and write to output 
  growseas_perc <- round((length(growseas_invalid)/nrow(data))*100,2)
  precip_perc   <- round((length(precip_invalid)/nrow(data))*100,2)
  PPFD_perc     <- round((length(PPFD_invalid)/nrow(data))*100,2)
  Tair_perc     <- round((length(Tair_invalid)/nrow(data))*100,2)
  ustar_perc    <- round((length(ustar_invalid)/nrow(data))*100,2)
  rH_perc       <- round((length(rH_invalid)/nrow(data))*100,2)
  
  cat(length(growseas_invalid)," data points (",growseas_perc,"%) excluded by growing season filter",fill=T,sep="")
  cat(length(precip_invalid)," data points (",precip_perc,"%) excluded by precipitation filter",fill=T,sep="")
  cat(length(PPFD_invalid)," data points (",PPFD_perc,"%) excluded by radiation filter",fill=T,sep="")
  cat(length(Tair_invalid)," data points (",Tair_perc,"%) excluded by air temperature filter",fill=T,sep="")
  cat(length(ustar_invalid)," data points (",ustar_perc,"%) excluded by friction velocity filter",fill=T,sep="")
  cat(length(rH_invalid)," data points (",rH_perc,"%) excluded by relative humidity filter",fill=T,sep="")
  
  invalid        <- unique(growseas_invalid,precip_invalid,PPFD_invalid,Tair_invalid,ustar_invalid,rH_invalid)
  valid[invalid] <- 0
  
  cat(length(invalid)," data points (",length(invalid)/nrow(data),"%) excluded in total",fill=T,sep="")

  # 5) set all other columns to NA
  data_filtered <- data.frame(data,valid)
  data_filtered[valid==0,!colnames(data) %in% c(precip,PPFD,Tair,ustar,VPD,GPP,day,month,year)] <- NA
  
  return(data_filtered)
}




#' GPP-based growing season filter
#' 
#' @description filters annual time series for growing season based on smoothed daily GPP time series
#' 
#' @param GPPd   daily GPP
#' @param tGPP   GPP threshold (fraction of 95% quantile) 
#' @param ws     window size used for smoothing
#' 
#' 
#' @details 



filter.growseas <- function(GPPd,tGPP,min_int,ws){  # ws = window size
  if(sum(is.na(GPPd)) < 0.5*length(GPPd)){
    growseas       <- rep(1,length(GPPd))
    GPP_threshold  <- quantile(GPPd,probs=0.95,na.rm=TRUE)*tGPP
    
    ## smooth GPP
    GPPd_smoothed      <- filter(GPPd,method="convolution",filter=rep(1/ws,ws))
    
    ## set edges to the mean of the original values
    wsd <- floor(ws/2)
    GPPd_smoothed[1:wsd] <- mean(GPPd[1:(2*wsd)],na.rm=T)
    GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] <- mean(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)],na.rm=T)
    
    # check for occurence of missing values and set them to mean of values surrounding them
    missing <- which(is.na(GPPd_smoothed))
    if (length(missing) > 0){
      if (length(missing) > 10){warning("Attention, there is a gap in GPPd of length n = ",length(missing))}
      replace_val <- mean(GPPd_smoothed[max(1,missing[1] - 4):min((missing[length(missing)] + 4),length(GPPd_smoothed))],na.rm=T)
      GPPd_smoothed[missing] <- replace_val
    }
    
    # filter daily GPP
    growseas[GPPd_smoothed < GPP_threshold] <- 0
    
    ## exclude short time periods
    intervals <- rle(growseas)
    short     <- which(intervals$lengths <= min_int)
    
    if (length(short) > 0){
      start <- numeric()
      end   <- numeric()
      for (i in 1:length(short)){
        start[i] <- sum(intervals$lengths[1:short[i]-1]) + 1
        end[i] <- start[i]+intervals$lengths[short[i]] - 1
        ## exclude
        val <- unique(growseas[start[i]:end[i]])
        if(length(val) > 1){
          stop("something went wrong when filtering growing season!")
          #       } else if (short[i] == 1){                          # only if first datapoints are concerned (use end instead of start!)
          #         if (val == 0 & growseas[end[1]+1] == 1){     
          #           growseas[start[1]:end[1]] <- growseas[start[1]:end[1]]+1
          #         } else if (val == 1 & growseas[end[1]+1] == 0){         
          #           growseas[start[1]:end[1]] <- growseas[start[1]:end[1]]-1
          #         } 
        } else {
          if (val == 0 & growseas[start[i]-1] == 1){              # second condition to avoid wrong fluctuations
            growseas[start[i]:end[i]] <- growseas[start[i]:end[i]]+1    # set to 1 (growing season)
          } else if (val == 1 & growseas[start[i]-1] == 0){
            growseas[start[i]:end[i]] <- growseas[start[i]:end[i]]-1    # set to 0 (no growing season)
          }
        }
      }
    }
    growseas <- as.integer(growseas)
  } else {
    warning("number of available GPPd data is less than half of total number of days per year. Filter is not applied!")
    growseas <- as.integer(rep(1,length(GPPd)))
  }
  growseas_invalid <- which(growseas==0)
  return(growseas_invalid)
}




###############################
### WUE metrics calculation ###-------------------------------------------------------------------
###############################


#' Water Use Efficiency metrics
#' 
#' @description Calculation of various water use efficiency metrics
#' 
#' @param data      Data.frame or matrix containing all required columns
#' @param GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' @param NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Tair      Air temperature (deg C)
#' @param constants 
#'
#' @details the following metrics are calculated:
#' 
#' @return a named vector with the following elements:
#' 
#' @references Beer_2009
#'             
#'   

### split up per year? or do this externally in aggregate function?!!! Is that working? 
### or would it be better to split it within the function
WUE.metrics <- function(data,GPP="GPP",NEE="NEE",LE="LE",VPD="VPD",Tair="Tair",
                        constants=bigleaf.constants()){
  
  GPP  <- check.columns(data,GPP)
  NEE  <- check.columns(data,NEE)
  LE   <- check.columns(data,LE)
  VPD  <- check.columns(data,VPD)
  Tair <- check.columns(data,Tair)
  
  
  ET  <- LEtoET(LE,Tair) * 86400          # kg H2O d-1
  GPP <- (GPP * 86400/1e06 * Cmol)*1000   # gC m-2 d-1
  NEE <- (NEE * 86400/1e06 * Cmol)*1000   # gC m-2 d-1
  
  
  WUE     <- median(GPP/ET,na.rm=T)
  WUE_NEE <- median(abs(NEE)/ET,na.rm=T)
  IWUE    <- median((GPP*VPD)/ET,na.rm=T)
  uWUE    <- median((GPP*sqrt(VPD))/ET,na.rm=T)
  
  return(c(WUE=WUE,WUE_NEE=WUE_NEE,IWUE=IWUE,uWUE=uWUE))
}



#' Estimation of g1
#' 
#' @description Estimation of the WUE metric "g1" from non-linear regression
#' 
#' @param data  a data.frame or matrix containing all required columns
#' @param GPP   gross primary productivity (umol m-2 s-1)
#' @param Gs    surface conductance (mol m-2 s-1)
#' @param VPD   vapor pressure deficit (kPa)
#' @param Ca    atmospheric CO2 concentration (ppm)
#' @param nmin  minimum number of data required to perform the fit
#' @param fitg0 logical, if TRUE g0 and g1 are fitted simultaneously
#' @param g0    ignored if fitg0 is FALSE, minimum stomatal conductance (mol m-2 s-1)
#' 
#' 
#' @return an nls model object
#' 
#' @references Medlyn 2011
#'             Knauer 2017
#' 
#' @export 

# basic tests on variable range?
estimate.g1 <- function(data,GPP="GPP",Gs="Gs",VPD="VPD",Ca="Ca",fitg0=FALSE,nmin=40,g0=0){
  GPP  <- check.columns(data,GPP)
  Gs   <- check.columns(data,Gs)
  VPD  <- check.columns(data,VPD)
  Ca   <- check.columns(data,Ca)
  
  nr_data <- sum(!is.na(GPP) & !is.na(Gs) & !is.na(VPD) & !is.na(Ca))
  

  if (nr_data < nmin){
    stop("number of data is less than 'nmin'. g1 is not fitted.")
  } else {
    if (fitg0){
      mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3))
    } else {
      mod <- nls(Gs ~ 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3))
    }
  } 
  
  return(mod)
    
}

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


################################################
### Boundary layer conductance formulations #### -------------------------------------------------------------------------
################################################

#' Boundary layer conductance according to Thom 1972
#' 
#' An empirical formulation for the canopy boundary layer conductance
#' based on a simple ustar dependency.
#' 
#' @param ustar     friction velocity (m s-1)
#' @param constants k - von-Karman constant (-)
#' 
#' @return a data.frame with the following columns:
#'  \item{Rb}{Boundary layer resistance (s m-1)}
#'  \item{Gb}{Boundary layer conductance (m s-1)}
#'  \item{kB}{kB-1 parameter (-)}
#' 
#' @references Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
#'             Quarterly Journal of the Royal Meteorological Society 98, 124-134.
#' 
#' 
#' @seealso \code{\link{Gb.Su}}, \code{\link{Gb.Choudhury}}
#' 
#' 
#' @export
Gb.Thom <- function(ustar,constants=bigleaf.constants()){
  Rb <- 6.2*ustar^-0.667
  Gb <- 1/Rb
  kB <- Rb*constants$k*ustar
  
  return(data.frame(Rb,Gb,kB))
}


#' Roughness Reynolds Number
#' 
#' @description 
#' 
#' 
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)
#' @param ustar     friction velocity (m s-1)
#' @param hs        roughness length of the soil (m)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K)
#'                  
#' @details The roughness Reynolds Number is calculated as (Massman 1999a):
#'          
#'          \deqn{Re = hs * ustar / v}
#'          
#'          where v is the kinematic viscosity of the air (m2 s-1), given by (Massman 1999b):
#'          
#'          \deqn{v = 1.327 * 10^-5(pressure0/pressure)(Tair/Tair0)^1.81}
#'        
#' @return a data.frame with the following columns:
#'         \item{v}{kinematic viscosity of air (m2 s-1)}
#'         \item{Re}{Roughness Reynolds Number (-)}
#' 
#' @references Massman, W.J., 1999a: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#' 
#'             Massman, W.J., 1999b: Molecular diffusivities of Hg vapor in air, 
#'             O2 and N2 near STP and the kinematic viscosity and thermal diffusivity
#'             of air near STP. Atmospheric Environment 33, 453-457.      
#'                 
#' 
#' @export
Reynolds.Number <- function(Tair,pressure,ustar,hs,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * 1000
  
  v  <- 1.327e-05*(constants$pressure0/pressure)*(Tair/constants$Tair0)^1.81
  Re <- hs*ustar/v
  
  return(data.frame("v"=v,"Re"=Re))
}
# hs = roughness height of the soil [m]
# v = kinematic viscosity of the air (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
# Tair in K!!




#' Boundary layer conductance according to Su et al. 2001
#' 
#' A physically based formulation for the canopy boundary layer conductance. 
#'
#' @param data      a matrix or data.frame containing all required variables
#' @param Tair      air temperature (degC)
#' @param presssure atmospheric pressure (Pa)
#' @param ustar     friction velocity (m s-1)
#' @param wind      wind speed (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param zh        canopy height (m)
#' @param Dl        leaf dimension (m)
#' @param N         number of leaf sides participating in heat exchange (1 or 2)
#' @param fc        fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI       one-sided leaf area index (-)
#' @param Cd        foliage drag coefficient (-)
#' @param hs        roughness height of the soil (m)
#' @param Pr        Prandtl number (-), defaults to 0.71
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K)
#' 
#' @return a data.frame with the following columns:
#'     \item{Rb}{Boundary layer resistance (s m-1)}
#'     \item{Gb}{Boundary layer conductance (m s-1)}
#'     \item{kB}{kB-1 parameter (-)}
#'     
#' @details The formulation is based on the kB-1 model developed by Massman 1999. 
#'          Su et al. 2001 derived the following approximation:
#'           
#'          \deqn{kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2
#'          
#'          If fc (fractional vegetation cover) is missing, it is estimated from LAI:
#' 
#'          \deqn(fc = 1 - exp(-LAI/2))
#'          
#'          The wind speed at the top of the canopy is calculated using function
#'          \code{\link{wind.profile}}.
#'          
#'          Ct is the heat transfer coefficient of the leaf (Massman 1999):
#'          
#'          \deqn{Ct = Pr^-2/3 Reh^-1/2 N}
#'          
#'          where Pr is the Prandtl number, and Reh is the Reynolds number for leaves:
#'          
#'          \deqn{Dl wind(zh) / v}
#'           
#'          kBs-1, the kB-1 value for bare soil surface, is calculated according 
#'          to Su et al. 2001:
#'          
#'          \deqn{kBs^-1 = 2.46(Re)^0.25 - ln(7.4)}
#' 
#' @references Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
#'             two models for estimation of the roughness height for heat transfer between
#'             the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
#' 
#'             Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#' 
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Choudhury}}
#' 
#' @export
Gb.Su <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",wind="wind",
                  H="H",zh,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,Pr=0.71,
                  constants=bigleaf.constants()){
        
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  wind     <- check.columns(data,wind)
  H        <- check.columns(data,H)
  
  if (is.null(fc)) {
    fc  <- (1-exp(-LAI/2)) 
  } 
  
  #wind_zh <- wind.profile(zh,data,) ## check first if that even works!!,check also if stab correctin is necessary
  
  v   <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"v"]
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"Re"]
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind / v
  Ct  <- 1*Pr^-0.6667*Reh^-0.5*N                       

  kB  <- (constants$k*Cd)/(4*Ct*ustar/wind)*fc^2 + kBs*(1 - fc)^2
  Rb  <- kB/(constants$k*ustar) 
  Gb  <- 1/Rb
  
  
  return(data.frame(Rb,Gb,kB))
}
# Ct = heat transfer coefficient of the leaf (Massman_1999 p.31)

# #' Boundary layer conductance according to McNaughton and van den Hurk 1995
# #' 
# #' @param ustar     friction velocity (m s-1)
# #' @param leafwidth leaf width (m)
# #' @param LAI       one-sided leaf area index
# #' @param constants k - von-Karman constant (-)
# #' 
# #' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Su}}
# 
# ## check: kB supposed to be so low??
# Gb.McNaughton <- function(ustar,leafwidth,LAI,constants=bigleaf.constants()){
#   Rb <- 130*(sqrt(leafwidth*ustar))/LAI - 1.7
#   
#   #kB <- k*(120/LAI*sqrt(leafwidth*ustar) - 2.5)
#   
#   kB  <- Rb*constants$k*ustar
#   Gb  <- 1/Rb
#   
#   
#   return(data.frame(Rb,Gb,kB))
# }



#' Boundary layer conductance according to Choudhury & Monteith 1988
#' 
#' @param wind      wind speed at sensor height (m s-1)
#' @param ustar     friction velocity (m s-1)
#' @param leafwidth leaf width (m)
#' @param LAI       one-sided leaf area index
#' @param zh        canopy height (m)
#' @param zr        sensor (instrument) height (m)
#' @param constants k - von-Karman constant (-)
#' 
#' @return a data frame with the following columns:
#'     \item{Rb}{Boundary layer resistance (s m-1)}
#'     \item{Gb}{Boundary layer conductance (m s-1)}
#'     \item{kB}{kB-1 parameter (-)}
#' 
#' @details Boundary layer conductance according to Choudhury & Monteith 1988 is
#'          given by:
#'          
#'          \deqn{Gb = LAI((2a/\alpha)*sqrt(u(h)/w)*(1-exp(-\alpha/2)))}
#'          
#'          where u(zh) is the wind speed at the canopy surface, approximated from
#'          measured wind speed at sensor height zr and a wind extinction coefficient \eqn{\alpha}:
#'          
#'          \deqn{u(zh) = u(zr) / (exp(\alpha(zr/zh -1)))}.
#'          
#'          \eqn{\alpha} is modeled as an empirical relation to LAI (McNaughton & van den Hurk 1995):
#'          
#'          \deqn{\alpha = 4.39 - 3.97*exp(-0.258*LAI)}
#'          
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             McNaughton, K. G., Van den Hurk, B.J.J.M., 1995: A 'Lagrangian' revision of
#'             the resistors in the two-layer model for calculating the energy budget of a
#'             plant canopy. Boundary-Layer Meteorology 74, 261-288.
#'             
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Su}}  
#'            
#' @export                                                                                                                                                                                                                                                                                    
Gb.Choudhury <- function(wind,ustar,leafwidth,LAI,zh,zr,constants=bigleaf.constants()){
  
  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)
  wind_zh <- wind / (exp(alpha*(zr/zh - 1)))
  Gb      <- LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  Rb      <- 1/Gb
  kB      <- Rb*constants$k*ustar
  
  return(data.frame(Rb,Gb,kB))
}


#' Wind speed at given levels above the canopy
#' 
#' @description wind speed at a given height above the canopy estimated from single-level
#'              measurements of wind speed at some distance above the canopy. 
#'              
#' @details The underlying assumption is the existence of an exponential wind profile
#'          above the height d + z0m (the height at which wind speed mathematically reaches zero
#'          according to the Monin-Obhukov similarity theory).
#'          In this case, the wind speed at a given height z is given by:
#'          
#'          \deqn{u(z) = (ustar/k) * (ln((z - d) / z0m) - \psi}
#'          
#'          The roughness parameters zero plane displacement height (d) and roughness length (z0m)
#'          can be estimated from \code{\link{roughness.parameters}}.
#'          
#' @param heights   vector with heights for which wind speed is to be 
#'                  calculated.
#' @param data      a matrix or data.frame containing all required variables
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)                                                                                  
#' @param ustar     friction velocity (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param zr        sensor (instrument) height (m)
#' @param d         displacement height (m)
#' @param z0m       roughness length for momentum (m)
#' @param constants k - von-Karman constant (-)
#'                                   
#' @note Note that this equation is only valid for z >= d + z0m. All values in \code{heights}
#'       smaller than d + z0m are set to NA.                                 
#'                                  
#' @return a data.frame with rows representing time and columns representing heights.     
#'                                            
#' @references                                                                                                                          

## deal with a) negative wind speeds
#            b) heights smaller than d + z0m
wind.profile <- function(heights,data,Tair="Tair",pressure="pressure",ustar="ustar",
                         wind="wind",H="H",zr,d,z0m,formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  formulation <- match.arg(formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  wind     <- check.columns(data,wind)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  
  psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d,
                                  formulation=formulation,constants)
  
  wind_heights <- data.frame(matrix(NA,ncol=length(heights),nrow=length(Tair)))
  colnames(wind_heights) <- paste0(heights,"m")
  for (z in heights){
    i <- which(heights == z)
    wind_heights[,i] <- (ustar / constants$k) * (log((z - d) / z0m) - psi_h)
  }
  
  return(wind_heights)
}


#' Roughness parameters
#' 
#' @description a simple approximation of the roughness parameters displacement height (d)
#'              and roughness length for momentum (z0m).
#'              
#' @param method    method used, either "canopy_height", or "wind_profile" \cr
#'                  NOTE: if method is "canopy_height", only the following three arguments
#'                  are used.               
#' @param frac_d    fraction of displacement height on canopy height (-)
#' @param frac_z0m  fraction of roughness length on canopy height (-)
#'                  The following arguments are only needed if method = "wind_profile"!
#' @param data      a matrix or data.frame containing all required variables
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)
#' @param wind      wind speed at height zr (m s-1)
#' @param ustar     friction velocity (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param zr        instrument height (m)
#' @param d         optional; displacement height (m)
#' @param z0m       optional; roughness length for momentum (m)
#' @param constants k - von-Karman constant (-)         
#' 
#' 
#' @details The two main roughness parameters, the displacement height (d)
#'          and the roughness length for momentum (z0m) can be estimated from simple
#'          empirical relationships with canopy height (zh). If \code{method} is
#'          \code{canopy_height}, the following formulas are used:  
#'          
#'          \deqn{d = frac_d * zh}
#'          
#'          \deqn{z0m = frac_z0m * zh}
#'          
#'          where frac_d defaults to 0.7 and frac_z0m to 0.1.
#'          
#'          If \code{method} is \code{wind_profile}, z0m is estimated by solving
#'          the windspeed profile for z0m:
#'          
#'                  
#' 
#' @return a data.frame with the following columns:
#'         \item{d}{displacement height (m)}
#'         \item{z0m}{roughness length for momentum (m)}
#'         


## estimation of d seems very unreliable! Evtl. exclude this case...
roughness.parameters <- function(method=c("canopy_height","wind_profile"),zh=NULL,frac_d=0.7,frac_z0m=0.1,data,
                                 Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",zr,d=NULL,z0m=NULL,
                                 formulation=c("Dyer_1970","Businger_1971"),constants=bigleaf.constants()){
  
  method      <- match.arg(method)
  formulation <- match.arg(formulation)
  
  if (method == "canopy_height"){
  
    if (is.null(zh)){
      stop("canopy height (zh) must be provided for method = 'canopy_height'!")
    }
    d   <- frac_d*zh
    z0m <- frac_z0m*zh
  
  } else if (method == "wind_profile"){
    
    Tair     <- check.columns(data,Tair)
    pressure <- check.columns(data,pressure)
    wind     <- check.columns(data,wind)
    ustar    <- check.columns(data,ustar)
    H        <- check.columns(data,H)
    
    if (is.null(d) & is.null(z0m)){
      
      d   <- frac_d*zh
      psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d,
                                      formulation=formulation,constants)
      
      z0m <- median((zr - d) * exp(-constants$k*wind / ustar - psi_h),na.rm=T)
    
    } else if (is.null(d) & !is.null(z0m)){
      
      d1    <-  frac_d*zh
      psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d=d1,
                                      formulation=formulation,constants)
      
      d <- median(zr - exp( (wind * constants$k) / ustar  - psi_h) * z0m,na.rm=T)
      
    } else if (!is.null(d) & is.null(z0m)){
      
      psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d=frac_d*zh,
                                      formulation=formulation,constants)
                                        
      z0m <- median((zr - d) * exp(-constants$k*wind / ustar - psi_h),na.rm=T)
    
    }
    
  }

  return(data.frame(d,z0m))
}


############################
### air pressure, density ## --------------------------------------------------------------------------
############################

#' Air density
#' 
#' @description Air density of moist air from air temperature and pressure
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (kPa)
#' @param constants Kelvin - conversion degC to Kelvin \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1)
#' 
#' @details Air density (\eqn{\rho}) is calculated as:
#' \deqn{\rho = pressure / Rd * Tair}
#' 
#' @return air density (kg m-3)
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 

air.density <- function(Tair,pressure,constants=bigleaf.constants()){
  Tair     <- Tair + constants$Kelvin
  pressure <- pressure*1000
  
  rho <- pressure / (constants$Rd * Tair) 
  
  return(rho)
}


#' Air pressure from hypsometric equation
#' 
#' @description An estimate of mean pressure at a given elevation as predicted by the
#'              hypsometric equation.
#'              
#' @param elev elevation a.s.l. (m)
#' @param Tair air temperature (degC)
#' @param q    specific humidity (kg kg-1); optional
#' @param constants Kelvin- conversion degC to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  g -  gravitational acceleration (m s-2) \cr
#' 
#' 
#' @details 
#' 
#' @note the hypsometric equation gives a rough estimate of the standard pressure
#'       at a given altitude. If specific humidity q is provided, humidity correction
#'       is applied and the virtual temperature instead of air temperature is used.
#'
#' @return air pressure (kPa)
#'                            
#' @references Stull B., 1988: An Introduction to Boundary Layer Meteorology.
#'             Kluwer Academic Publishers, Dordrecht, Netherlands.
#'             
pressure.from.elevation <- function(elev,Tair,q=NULL,constants=bigleaf.constants()){
  if(is.null(q)){
    Temp <- Tair + constants$Kelvin
  } else {
    Temp <- Temp.virtual(Tair,q,constants)  + constants$Kelvin
  }
  
  pressure <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Temp))
  pressure <- pressure/1000
  
  return(pressure)
} 



















###########################################
### Stability parameters and correction ### ---------------------------------------------------
###########################################

#' Monin-Obukhov length
#' 
#' @description 
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (kPa)
#' @param ustar     friction velocity (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  k - von Karman constant (-) \cr
#'                  g - gravitational acceleration (m s-2)
#' 
#' @description The Monin-Obukhov length (L) is given by:
#' 
#'              \deqn{L = - (\rho * cp * ustar^3 * Tair) / (k * g * H)}
#' 
#' @return Monin-Obukhov length (m)
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
#' @export
MoninObukhov.length <- function(Tair,pressure,ustar,H,constants=bigleaf.constants()){
  
  rho  <- air.density(Tair,pressure,constants=bigleaf.constants())
  Tair <- Tair + constants$Kelvin
  MOL  <- (-rho*constants$cp*ustar^3*Tair) / (constants$k*constants$g*H)
  
  return(MOL)
}


#' stability parameter 'zeta'
#' 
#' @param Tair      air temperature (degC)
#' @param pressure  air pressure (kPa)
#' @param ustar     friction velocity (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param zr        instrument height (m)
#' @param disp      displacement height (m)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  k - von Karman constant (-) \cr
#'                  g - gravitational acceleration (m s-2)
#' 
#' @details The stability parameter \eqn{\zeta} is given by:
#' 
#'          \deqn{\zeta = (zr - d) / L}
#'          
#'          where L is the Monin-Obukhov length, calculated from the function
#'          \code{\link{MoninObukhov.length}}. The displacement height d can 
#'          be estimated from the function \code{\link{roughness.parameters}}.
#'          
#' @return the stability parameter \eqn{\zeta}
#' 
#' @export           
stability.parameter <- function(Tair,pressure,ustar,H,zr,d,constants=bigleaf.constants()){
  
  MOL  <- MoninObukhov.length(Tair,pressure,ustar,H,constants)
  zeta <- (zr - d) / MOL
  
  return(zeta)
  
}


#' stability correction functions for sensible heat
#' 
#' @description dimensionless stability functions needed to correct deviations
#'              from the exponential wind profile under non-neutral conditions.
#'              
#' @param zeta         stability parameter \eqn{\zeta}
#' @param formulation  formulation for the stability function. Either "Dyer_1970", or "Businger_1971"
#'  
#' @details The functions depend on the value of the stability parameter \eqn{\zeta},
#'          which can be calculated from the function \code{\link{stability.parameter}}.
#'          The general form of the stability functions is:
#'          
#'          \deqn{-x * zeta} 
#'          
#'          for stable atmospheric conditions (\eqn{\zeta} >= 0), and
#'          
#'          \deqn{2 * log( (1 + y) / 2) }
#'          
#'          for unstable atmospheric conditions (\eqn{\zeta} < 0).
#'          
#'          The different formulations differ in their value of x and y.
#'   
#' @return psi_h the value of the stability function for sensible heat (-)  
#' 
#' @references Dyer, A.J., 1974: A review of flux-profile relationships. 
#'             Boundary-Layer Meteorology 7, 363-372.
#'             
#'             Dyer, A. J., Hicks, B.B., 1970: Flux-Gradient relationships in the
#'             constant flux layer. Quart. J. R. Meteorol. Soc. 96, 715-721.
#'             
#'             Businger, J.A., Wyngaard, J. C., Izumi, I., Bradley, E. F., 1971:
#'             Flux-Profile relationships in the atmospheric surface layer. 
#'             J. Atmospheric Sci. 28, 181-189.
#' 
#'             Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#'   

stability.correction.H <- function(Tair,pressure,ustar,H,zr,d,formulation=c("Dyer_1970","Businger_1971"),
                                   constants=bigleaf.constants()){
  
  formulation  <- match.arg(formulation)
  
  zeta <- stability.parameter(Tair,pressure,ustar,H,zr,d,constants)
  
  psi_h <- numeric()
  
  if (formulation == "Businger_1971"){
    x <- -7.8
    y <- 0.95 * ( 1 - 11.6 * zeta)^0.5
  } else if (formulation == "Dyer_1970"){
    x <- -5
    y <- (1 - 16 * zeta)^0.5
  }

  # stable
  psi_h[zeta >= 0] <- -x * zeta[zeta >= 0]
  # unstable
  psi_h[zeta < 0] <- 2 * log(( 1 + y[zeta < 0] ) / 2)  
  
  return(psi_h)
} 

### Dyer_1970; Dyer_1974 (similar:Paulson 1970)
### Businger_1971: from Foken_2008; p.65 (suggested by Businger_1971), in the form of Högstrom 1988





#' Aerodynamic conductance
#' 
#' @description Aerodynamic conductance, including options for the boundary layer conductance
#'              formulation and stability correction functions.
#' 
#' @param Tair              air temperature (degC)
#' @param pressure          air pressure (Pa)
#' @param wind              wind speed (m s-1)
#' @param ustar             friction velocity (m s-1)
#' @param H                 sensible heat flux (W m-2)
#' @param zr                measurement (=reference) height (m)
#' @param zh                canopy height (m)
#' @param d                 zero-plane displacement height (m)
#' @param z0m               roughness length for momentum (m)
#' @param Dl                characteristic leaf dimension (m)
#' @param N                 number of leaf sides participating in heat exchange (1 or 2)
#' @param LAI               one-sided leaf area index (m2 m-2)
#' @param fc                fractional vegetation cover (-)
#' @param Cd                foliage drag coefficient (-)
#' @param hs                roughness length of bare soil (m)
#' @param stab_correction   stability correction function
#' @param Rb_model          boundary layer resistance formulation
#' @param kB                kB-1 value; only used if \code{Rb_model = "constant_kB"}
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
#'  where z = reference height, d the zero-plane displacement height, and L the Monin-Obukhov length, 
#'  calculated with \code{\link{MoninObukhov.length}}
#'  Stability correction functions calculate the correction based on the formulations by Businger 1971 
#'  Dyer 1970, based on zeta.
#'  If "none", atmospheric stability is neglected, and Monin-Obukhov length, as well as the 
#'  stability parameter zeta are not calculated.
#' 
#'  The argument \code{Rb_model} determines the canopy boundary layer resistance ("excess resistance") model that is used. "Thom_1972" is an 
#'  empirical formulation based on the friction velocity (ustar) (Thom 1972):
#'  
#'    \deqn{Rb = 0.667u*^0.667}
#'    
#'  The model by Mcnaughton and van den Hurk 1995 ("McNaughton_1995"), calculates Rb
#'  based on leaf width, LAI and ustar (Note that the original formulation leaf width
#'  instead of the characteristic leaf dimension (Dl) is taken):
#'   
#'     formula
#'     
#'  The option "Su_2001" calculates Rb based on the physically-based Rb model by Su et tal. 2001,
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
#' @references 
#'  
#' @export
#'    


# TO DO:
# - find out what "." means in front of a function!)
# - think about how to deal with uncertainty!
#- think about testing for missing arguments in the subfunctions separately, but most likely doesn't make sense
Ga <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",zr,zh,d,z0m,
               frac_d=0.7,frac_z0m=0.1,Dl,N,LAI,fc=NULL,Cd=0.2,hs=0.01,
               stab_correction=c("Dyer_1970","Businger_1971","none"),
               method_roughness=c("canopy_height","wind_profile"),
               Rb_model=c("Thom_1972","Choudhury_1988","Su_2001","constant_kB"),
               kB=NULL,constants=bigleaf.constants()) 
      {
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  wind     <- check.columns(data,wind)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  
  Rb_model         <- match.arg(Rb_model)
  stab_correction  <- match.arg(stab_correction)
  method_roughness <- match.arg(method_roughness)
  
  
  if (Rb_model != "constant_kB"){
    
    if (Rb_model == "Thom_1972"){
      
      Gb_mod <- Gb.Thom(ustar,constants)
      
    } else if (Rb_model == "Choudhury_1988"){
      
      Gb_mod <- Gb.Choudhury(wind,ustar,Dl,LAI,zh,zr,constants)
      
    } else if (Rb_model == "Su_2001"){
      
      Gb_mod <- Gb.Su(ustar,wind,pressure,Tair,Dl,N,fc=fc,LAI,Cd=Cd,hs=hs,constants)
      
    }
 
      kB <- Gb_mod[,"kB"]
      Rb <- Gb_mod[,"Rb"]
      Gb <- Gb_mod[,"Gb"]
  
  } else {
    
    if(is.null(kB)){
      stop("value of kB has to be specified if Rb_model is set to 'constant_kB'!")
    } else {
      Rb <- kB / (constants$k * ustar) 
      Gb <- 1 / Rb
    }
  
  }
  
  Ra_m <- wind/ustar^2
  
  
  if (stab_correction != "none"){
      
    rp <- roughness.parameters(method=method_roughness,zh,frac_d,frac_z0m,data,zr=zr,
                               d=d,z0m=z0m,formulation=stab_correction,constants)
    d   <- rp[,"d"]
    z0m <- rp[,"z0m"]
    
    
    zeta  <-  stability.parameter(Tair,pressure,ustar,H,zr,d,constants)
    
    
    if (stab_correction == "Dyer_1970"){
      
      psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d,formulation="Dyer_1970")
      
    } else if (stab_correction == "Businger_1971"){
      
      psi_h <- stability.correction.H(Tair,pressure,ustar,H,zr,d,formulation="Businger_1971")
      
    }  
    
    Ra_m  <- pmax((log((zr - d)/z0m) - psi_h),1e-10)/(constants$k*ustar)
  
  
  } else {
    
    zeta = psi_h <- rep(NA,length=length(Ra_m))
    
  }
  
  
  Ga_m <- 1/Ra_m
  Ra_h <- Ra_m + Rb
  Ga_h <- 1/Ra_h
   

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Rb=Rb,Gb=Gb,kB=kB,zeta,psi_h))
  
}
  

### alternative approaches for estimating d and z0m
#     diff <- numeric(); ustar1 <- ustar; z0m <- numeric()
#     for (i in 1:95){
#       ustar1[ustar < quantile(ustar,i/100,na.rm=T)] <- NA
#       z0m[i]          <- summary(nls(Ra_m ~ (log((zr - disp)/z0m)/(constants$k*ustar1)),
#                                 start=c(z0m=z0m_start),na.action=na.omit,
#                                 data=environment()))$par[1]
#       
#       diff[i] <- sum((Ra_m - (log((zr - disp)/z0m[i])/(constants$k*ustar1)))^2,na.rm=T)
#     }

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



#' Surface conductance for water vapor
#' 
#' @description Calculates surface conductance for water vapor from the inverted Penman-Monteith
#'              equation.
#' 
#' @param data      data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         optional, Ground heat flux (W m-2)
#' @param S         optional, Sum of all storage fluxes (W m-2)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance (m s-1)
#' @param missing.G.as.NA  if TRUE, missing G are treated as NA,otherwise set to 0. 
#' @param missing.S.as.NA  if TRUE, missing S are treated as NA,otherwise set to 0. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin
#' 

#' 
#' @details Surface conductance is calculated from the inverted Penman-Monteith equation:
#' 
#'  formula
#'  
#'  Note that Gs > Gc (canopy conductance) in time periods when a significant fraction of 
#'  ET comes from interception or soil evaporation. 
#'  Available energy A is given by (Rn - G - S). If G and/or S are not provided, then A = Rn.
#'  
#'  By default, any missing data in G and S are set to 0. If arguments \code{missing.S.as.NA} or \code{missing.S.as.NA} are 
#'  set to TRUE, Gs will give NA for these timesteps.
#'  
#'  If pressure is not available, it can be approximated by elevation using the 
#'  function \code{\link{pressure.from.elevation}}
#'  
#' 
#' @return a dataframe with two columns: 
#'  \item{Gs_ms}{Surface conductance in ms-1}
#'  \item{Gs_mol}{Surface conductance in mol m-2 s-1}
#' 
#' @export
GsPM <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,VPD="VPD",LE="LE",Ga="Ga",
                 constants=bigleaf.constants(),missing.G.as.NA=FALSE,missing.S.as.NA=FALSE){
   
   Tair     <- check.columns(data,Tair)
   pressure <- check.columns(data,pressure)
   Rn       <- check.columns(data,Rn)
   VPD      <- check.columns(data,VPD)
   LE       <- check.columns(data,LE)
   Ga       <- check.columns(data,Ga)
  
   if(!is.null(G)){
     G <- check.columns(data,G)
     if (!missing.G.as.NA){G[is.na(G)] <- 0}
   } else {
     warning("ground heat flux G is not provided and set to 0.")
     G <- rep(0,nrow(data))
   }
   
   if(!is.null(S)){
     S <- check.columns(data,S)
     if(!missing.S.as.NA){S[is.na(S)] <- 0 }
   } else {
     warning("Energy storage fluxes S are not provided and set to 0.")
     S <- rep(0,nrow(data))
   }
   
  
   delta <- Esat(Tair)[,"delta"]
   gamma <- psychrometric.constant(Tair,pressure)
   rho   <- air.density(Tair,pressure)
 
   Gs_ms  <- ( LE * Ga * gamma ) / ( delta * (Rn-G-S) + rho * constants$cp * Ga * VPD - LE * ( delta + gamma ) )
   Gs_mol <- mstomol(Gs_ms,Tair,pressure) 

   return(data.frame(Gs_ms,Gs_mol))

}


#' Big-leaf surface conditions
#' 
#' @description Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' @param data      data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param LE        Latent heat flux (W m-2)
#' @param H         Sensible heat flux (W m-2)
#' @param Ca        Atmospheric CO2 concentration (mol mol-1)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance (m s-1)          
#' @param NEE       Net ecosystem exchange (umol m-2 s-1)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#' @details 
#' 
#' @return a data.frame with the following columns:
#'         \item{Tsurf}{Surface temperature (deg C)} \cr
#'         \item{esat_surf}{Saturation vapor pressure at the surface (kPa)} \cr
#'         \item{esurf}{vapor pressure at the surface (kPa)} \cr
#'         \item{VPD_surf}{vapor pressure deficit at the surface (kPa)} \cr
#'         \item{qsurf}{specific humidity at the surface} \cr
#'         \item{Ca_surf}{CO2 concentration at the surface}             
#'         
#'                                       
#' @export 
#' 

# - check Ga
# - recalculate Ga for Ca surface
# e too high
bigleaf.surface <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",Ca="Ca",
                            VPD="VPD",Ga="Ga",NEE="NEE",constants=bigleaf.constants()){

  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  LE       <- check.columns(data,LE)
  H        <- check.columns(data,H)
  Ca       <- check.columns(data,Ca)
  VPD      <- check.columns(data,VPD)
  Ga       <- check.columns(data,Ga)
  NEE      <- check.columns(data,NEE)
  
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
  qsurf     <- VPD.to.q(VPD,Tair,pressure,constants)

  # 3) CO2 concentration
  Ca_surf <- Ca + NEE/Ga

  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,Ca_surf))
}



## check: either Smith 1998 or Daudet are wrong? typo? Daudet is identical to Jarvis1986 at least. 
#' Decoupling coefficient
#' 
#' @description The canopy-atmosphere decoupling coefficient 'Omega'. 
#' 
#' @param data       data.frame or matrix containing all required input variables
#' @param Tair       air temperature (deg C)
#' @param pressure   air pressure (kPa)
#' @param Ga         aerodynamic conductance (m s-1)
#' @param Gs         surface conductance (m s-1)
#' @param approach   approach used to calculate omega. One of "JarvisMcNaughton_1986",
#'                   "Smith_1998","Daudet_1999" or "Martin_1989".
#' @param constants  Kelvin - conversion degree Celsius to Kelvin \cr
#'                   sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                   cp - specific heat of air for constant pressure (J K-1 kg-1) 
#' @param LAI        leaf area index (m2 m-2), only used if approach = "Martin_1989"
#' 
#' @details The decoupling coefficient Omega ranges from 0 to 1 and quantifies the
#'          linkage of the conditions (foremost humidity) at the canopy surface
#'          to the ambient air. Values close to 0 indicate well coupled conditions
#'          characterized by high physiological (i.e. stomatal) control on transpiration
#'          and similar conditions at the canopy surface compared to the atmosphere above
#'          the canopy. Values close to 1 indicate the opposite, i.e. a low physiological 
#'          control on transpiration (Jarvis & McNaughton 1986). \cr
#'          The "JarvisMcNaughton_1986" approach (default option) is the original
#'          decoupling coefficient, given by:
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1}{\epsilon + 1 + \frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 1) / ( \epsilon + 1 + Ga/Gc)}
#'          
#'          where \eqn{\epsilon = \frac{s}{\gamma}}{\epsilon = s/\gamma} is a dimensionless coefficient
#'          with s being the slope of the saturation vapor pressure curve (Pa K-1), and \eqn{\gamma} the 
#'          psychrometric constant (Pa K-1).
#'          
#'          The approach "Smith_1998" by Smith et al. 1998 has been proposed for 
#'          application on hypostomatous vegetation: 
#'          
#'          \deqn{\Omega = \frac{\epsilon + 2}{\epsilon + 2 + \frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 2) / (\epsilon + 2 + Ga/Gc)}
#'          
#'          Alternatively, "Daudet_1999" by Daudet et al. 1999 propose the following
#'          formulation for hypostomatous vegetation: 
#'          
#'          \deqn{\Omega = \frac{\epsilon + 2}{\epsilon + 2 + 2\frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 2) / (\epsilon + 2 + 2Ga/Gc)}
#'          
#'          The approach "Martin_1989" by Martin 1989 additionally takes radiative coupling
#'          into account:
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1 + \frac{Gr}{Ga}}{\epsilon + (1 + \frac{Ga}{Gs}) (1 + \frac{Gr}{Ga})}}{%
#'          \Omega = (\epsilon + 1 + Gr/Ga) / (\epsilon + (1 + Ga/Gs) (1 + Gr/Ga))}
#' 
#' @return the decoupling coefficient Omega (-)
#' 
#' @references Goldberg V., Bernhofer C., 2008: Testing different decoupling coefficients
#'             with measurements and models of contrasting canopies and soil water conditions.
#'             Annales Geophysicae 26, 1977-1992.
#' 
#'             Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49. 
#'             
#'             Smith D. M., Jarvis P. G., Odongo J. C. W, 1998: Management of
#'             windbreaks in the Sahel: the strategic implications of tree water
#'             use. Agroforestry Systems 40, 83-96.
#'             
#'             Daudet F. A., Le Roux X., Sinoquet H., and Adam A., 1999:
#'             Wind speed and leaf boundary conductance variation within
#'             tree crown, Consequences on leaf-atmosphere coupling and tree
#'             functions, Agricultural and Forest Meteorology 97, 171-185.
#'             
#'             Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#' 
#' @export
decoupling <- function(data,Tair="Tair",pressure="pressure",Ga="Ga",Gs="Gs",
                       approach=c("JarvisMcNaughton_1986","Smith_1998","Daudet_1999","Martin_1989"),
                       constants=bigleaf.constants(),LAI=NULL){

  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  Ga       <- check.columns(data,Ga)
  Gs       <- check.columns(data,Gs)

  delta   <- Esat(Tair)[,"delta"]
  gamma   <- psychrometric.constant(Tair,pressure,constants=constants)
  epsilon <- delta/gamma
   
  approach <- match.arg(approach)
  
  if (approach == "JarvisMcNaughton_1986"){
    
    omega <- (epsilon + 1) / (epsilon + 1 + Ga/Gs)
    
  } else if (approach == "Smith_1998"){
    
    omega <- (epsilon + 2) / (epsilon + 2 + Ga/Gs)  
    
  } else if (approach == "Daudet_1999"){
    
    omega <- (epsilon + 2) / (epsilon + 2 + 2*Ga/Gs) 
    
  } else if (approach == "Martin_1989") {
  
    if (is.null(LAI)){
      stop("LAI is not provided!")
    } else {
      Gr    <- Gr.longwave(Tair,LAI,constants=constants)
      omega <- (epsilon + 1 + Gr/Ga) / (epsilon + 1 + Ga/Gs + Gr/Gs + Gr/Ga)
    }
  }
  
  return(omega)

}


#' longwave radiative transfer conductance of the canopy
#' 
#' @param Tair      air temperature (deg C)
#' @param LAI       leaf area index (m2 m-2)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1)
#'                  
#' @details the following formula is used (Martin, 1989):
#' 
#'          \deqn{Gr = 4 \sigma Tair^3 LAI / cp}                             
#'                                       
#' @return Gr longwave radiative transfer conductance of the canopy (m s-1)                 
#'                  
#' @references Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#'             
#' @export             
Gr.longwave <- function(Tair,LAI,constants=bigleaf.constants()){
  
  Tair <- Tair + constants$Kelvin
    
  Gr <- 4 * constants$sigma * Tair^3 * LAI / constants$cp
  
  return(Gr)
}




#' virtual temperature
#' 
#' @description virtual temperature
#' 
#' @param Tair      air temperature (deg C)
#' @param q         specific humidity (kg kg-1)
#' @param constants Kelvin - conversion degree Celsius to Kelvin
#' 
#' @details the virtual temperature is given by:
#'  \deqn{Tv = (Tair+273.15)(1 + 0.61*q)}
#' 
#' @note mixing ratio is approximated by specific humidity.
#' 
#' @return virtual temperature (deg C)
#' 
#' @references Foken T., 2008: Micrometeorology. Springer, Berlin, Germany.
#' 
#' @export
Temp.virtual <- function(Tair,q,constants=bigleaf.constants()){
   Tair <- Tair + 273.15
   Tv   <- Tair*(1+0.61*q)   # mixing ratio is approximated by specific humidity
   Tv   <- Tv - constants$Kelvin
   return(Tv)
}


#' Saturation vapor pressure (Esat) and slope of the Esat curve
#' 
#' @description Calculates Saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' @param Tair     air temperature (deg C)
#' @param formula  formula to be used. Either "Alduchov_1996" or "Sonntag_1990"
#' 
#' @details Esat (kPa) is calculated based on the Magnus equation:
#' 
#'  \deqn{Esat = a * exp((b * Tair) / (c + Tair)) / 1000}
#'  
#'  where the coefficients a, b, c take different values depending on the formula used.
#'  The slope of the Esat curve is calculated as the first derivative of the function:
#'  
#'  formula
#' 
#' @return A dataframe with the following columns: 
#'    \item{Esat}{Saturation vapor pressure (kPa)}
#'    \item{delta}{Slope of the saturation vapor pressure curve (kPa K-1)}
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
  Esat <- Esat/1000
  
  delta <- eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  delta <- delta/1000
  
  
  return(data.frame(Esat,delta))
}



#' Psychrometric constant
#' 
#' @description Calculates the psychrometric 'constant'.
#' 
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' @details The psychrometric constant (\eqn{\gamma}) is given as:
#' 
#'  \deqn{\gamma = cp * pressure / (eps * \lambda)},
#'  
#'  where \eqn{\lambda} is the latent heaf of vaporization (J kg-1), 
#'  as calculated with \code{\link{LE.vaporization}}.
#'  
#' @return the psychrometric constant (kPa K-1)
#'  
#' @references Bonan_2008 or Monteith_2010
#' 
psychrometric.constant <- function(Tair,pressure,constants=bigleaf.constants()){
  
  lambda <- LE.vaporization(Tair)
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
#' @return \eqn{\lambda} latent heat of vaporization (J kg-1) 
#' 
#' @references Stull, B., 1988: An Introduction to Boundary Layer Meteorology (p.641)
#'             Kluwer Academic Publishers, Dordrecht, Netherlands
#'             
#' @export
LE.vaporization <- function(Tair) {
  k1 <- 2.501
  k2 <- 0.00237
  lambda <- ( k1 - k2 * Tair ) * 1e+06
  
  return(lambda)
}



#######################
## Unit Conversions ### ---------------------------------------------------------------------------
#######################

#' Conversion between latent heat flux and evapotranspiration
#' 
#' @description converts evpaorative water flux from mass (ET=evapotranspiration)
#'              to energy (LE=latent heat) units, or vice versa.
#'              
#' @aliases LEtoET ETtoLE
#' 
#' @param LE   Latent heat flux (W m-2)
#' @param ET   Evapotranspiration (kg m-2 s-1)
#' @param Tair Air temperature (deg C)
#' 
#' @export
LEtoET <- function(LE,Tair){

  Tair   <- Tair
  lambda <- LE.vaporization(Tair)
  ET     <- LE/lambda
   
  return(ET)
}

ETtoLE <- function(ET,Tair){
  
  Tair   <- Tair
  lambda <- LE.vaporization(Tair)
  LE     <- ET*lambda
  
  return(LE)
}


#' Conversion between conductance units
#' 
#' @description Converts canopy/surface conductance from ms-1
#'              to mol m-2 s-1, or vice versa.
#' 
#' @aliases mstomol moltoms
#' 
#' @param Gs_ms       Canopy/Surface conductance (ms-1)
#' @param Gs_mol      Canopy/Surface conductance (mol m-2 s-1)
#' @param Tair        Air temperature (deg C)
#' @param pressure    Air pressure (kPa)
#' @param constants   Kelvin - conversion degree Celsius to Kelvin \cr
#'                    Rgas - universal gas constant (J mol-1 K-1)
#' 
#' @family conductance conversion
#' 
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology.
#'             2nd Edition., 2nd Edn. Cambridge University Press, Cambridge. 428 p --> replace with third edition
#'             
#' @export
mstomol <- function(Gc_ms,Tair,pressure,constants=bigleaf.constants()){
  Tair   <- Tair + constants$Kelvin
  Gc_mol <- Gc_ms * pressure*1000 / (constants$Rgas * Tair)
  
  return(Gc_mol)
}


#' @rdname mstomol
#' @family conductance conversion
moltoms <- function(Gc_mol,Tair,pressure,constants=bigleaf.constants()){
  Tair  <- Tair + constants$Kelvin
  Gc_ms <- Gc_mol * (constants$Rgas * Tair) / (pressure*1000)
  
  return(Gc_ms)
}



#' Conversions between humidity measures
#' 
#' @description Conversion between vapor pressure, vapor pressure deficit (VPD), specific humidity, and relative humidity.
#' 
#' @param e         vapor pressure (kPa)
#' @param p         air pressure (kPa)
#' @param q         specific humidity (kg kg-1)
#' @param VPD       vapor pressure deficit (kPa)
#' @param rH        relative humidity (-)
#' @param constants eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @family humidity conversion
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#' 
#' @export
e.to.q <- function(e,pressure,constants=bigleaf.constants()){
  q <- constants$eps * e / (pressure - (1-constants$eps) * e) 
  return(q)
}


#' @rdname e.to.q
#' @family humidity conversion
#' @export
q.to.e <- function(q,pressure,constants=bigleaf.constants()){
  e <- q * pressure / ((1-constants$eps) * q + constants$eps)
  return(e)
}


#' @rdname e.to.q
#' @family humidity conversion
#' @export
e.to.VPD <- function(e,pressure,Tair){
  esat <- Esat(Tair)[,"Esat"]
  VPD  <- esat - e 
  return(VPD)
}


#' @rdname e.to.q
#' @family humidity conversion
#' @export
VPD.to.e <- function(VPD,pressure,Tair){
  esat <- Esat(Tair)[,"Esat"]
  e    <- esat - VPD
  return(e)
}


#' @rdname e.to.q
#' @family humidity conversion
#' @export
rH.to.VPD <- function(rH,Tair){
  esat <- Esat(Tair)[,"Esat"]
  VPD  <- esat - rH*esat
  return(VPD)
} 


#' @rdname e.to.q
#' @family humidity conversion
#' @export
VPD.to.rH <- function(VPD,Tair){
  esat <- Esat(Tair)[,"Esat"]
  rH   <- 1 - VPD/esat
  return(rH)
} 


#' @rdname e.to.q
#' @family humidity conversion
#' @export
q.to.VPD <- function(q,Tair,pressure,constants=bigleaf.constants()){
  esat <- Esat(Tair)[,"Esat"]
  e    <- q.to.e(q,pressure,constants)
  VPD  <- esat - e
  return(VPD)
} 


#' @rdname e.to.q
#' @family humidity conversion
#' @export
VPD.to.q <- function(VPD,Tair,pressure,constants=bigleaf.constants()){
  esat <- Esat(Tair)[,"Esat"]
  e    <- esat - VPD
  q    <- e.to.q(e,pressure,constants)
  return(q)
} 





#' radiation unit conversion ?
#'
#'
#'
moltoWatt <- function(PPFD,conv.factor){
  PAR <- PPFD * conv.factor
  return(PAR)
}

Watttomol <- function(PAR,conv.factor){
  PPFD <- PAR
  return(PPFD)
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
# 
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

###########################
#### Evapotranspiration ###-------------------------------------------------------------------------------
###########################

#' Potential evapotranspiration from the Priestley-Taylor equation
#' 
#' @details 
#' 
#' @param data      data.frame or matrix containing all required columns
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)
#' @param Rn        net radiation (W m-2)
#' @param G         ground heat flux (W m-2)
#' @param S         sum of all storage fluxes (W m-2)
#' @param alpha     Priestley-Taylor coefficient (-)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @details potential Evapotranspiration is calculated according to Priestley & Taylor, 1972:
#' 
#'          \deqn{(\alpha * \delta * (Rn - G - S)) / (\delta + \gamma)}
#'
#' @return ET_pot potential evapotranspiration (kg m-2 s-1) 
#'   
#' @references Priestley C.H.B., Taylor R.J., 1972: On the assessment of surface heat flux
#'             and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.           
#' 
#' @export
ET.pot <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,alpha=1.26,
                   missing.G.as.NA=F,missing.S.as.NA=F,constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  Rn       <- check.columns(data,Rn)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,nrow(data))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,nrow(data))
  }
  
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  delta  <- Esat(Tair)[,"delta"]
  
  LE_pot <- (alpha * delta * (Rn - G - S)) / (delta + gamma)
  ET_pot <- LEtoET(LE_pot,Tair)
  
  return(ET_pot)
}




#' Imposed and Equilibrium Evapotranspiration
#' 
#' @details Evapotranspiration (ET) split up in imposed ET and equilibrium ET.
#' 
#' @param data      data.frame or matrix containing all required input variables
#' @param Tair      air temperature (deg C)
#' @param pressure  air pressure (kPa)
#' @param Gs        surface conductance (m s-1)
#' @param VPD       air vapor pressure deficit (kPa)
#' @param Rn        net radiation
#' @param constants 
#' 
#' @return 
#' 
#' @references Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49.
#'             
#'             Monteith J.L., Unsworth M.H., 2008: Principles of environmental physics.
#'             3rd edition. Academic Press, London. 

## should Rn be available energy????
## double check this function again!
ET.components <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",
                          Rn="Rn",constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  Rn       <- check.columns(data,Rn)
  VPD      <- check.columns(data,VPD)
  
  rho    <- air.density(Tair,pressure)
  lambda <- LE.vaporization(Tair)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  delta  <- Esat(Tair)[,"delta"]
  
  ET_imp <- (rho * constants$cp * Gs * VPD) / gamma
  ET_equ <- (delta * Rn) / (gamma + delta)
  
  return(data.frame(ET_imp,ET_equ))
}



## ET from Penman-Monteith and Priestley-Taylor



  

########################
#### Energy balance ####-------------------------------------------------------------------------------
########################

#' Biochemical energy
#' 
#' @description Radiant energy absorbed in photosynthesis or heat release by respiration calculated
#'              from net ecosystem exchange of CO2 (NEE).  
#' 
#' @param alpha   energy taken up/released by photosynthesis/respiration (J umol-1)
#' @param NEE     net ecosystem exchange (umol CO2 m-2 s-1)
#' 
#' @details the following sign convention is employed: NEE is negative when carbon is taken up by the ecosystem.
#'          Positive values of the resulting biochemical energy mean that energy is taken up by the ecosystem, 
#'          negative ones that heat is released.
#'          The value of alpha is taken from Nobel 1974 (see Meyers & Hollinger 2004), but other values
#'          have been used (e.g. Blanken et al., 1997)
#' 
#' @return biochemical energy Sp (W m-2)
#' 
#' @references Meyers, T.P., Hollinger, S.E. 2004: An assessment of storage terms in the surface energy
#'             balance of maize and soybean. Agricultural and Forest Meteorology 125, 105-115.
#'             
#'             Nobel, P.S., 1974: Introduction to Biophysical Plant Physiology.
#'             Freeman, New York.
#'             
#'             Blanken, P.D. et al., 1997: Energy balance and canopy conductance of a boreal aspen
#'             forest: Partitioning overstory and understory components. 
#'             Journal of Geophysical Research 102, 28915-28927. 
#'            
#' @export 
biochemical.energy <- function(alpha=0.422,NEE){
  Sp <- -alpha*NEE
  return(Sp)
}



#' Energy balance closure
#' 
#' @description Calculates the degree of the energy balance non-closure based on the ratio
#'              of two sums (energy balance ratio), or ordinary least squares (OLS).
#' 
#' @param data data.frame or matrix containing all required columns
#' @param Rn   net radiation (W m-2)
#' @param G    ground heat flux (W m-2)
#' @param S    sum of all storage fluxes (W m-2)
#' @param LE   latent heat flux (W m-2)
#' @param H    sensible heat flux (W m-2)
#' 
#' 
#' @details The energy balance ratio (EBR) is calculated as:
#'          
#'          \deqn{EBR = sum(LE + H)/sum(Rn - G - S)}
#'          
#'          the sum is taken for all timesteps with complete observations (i.e. where
#'          all energy balance terms are available).
#' 
#' @return a named vector containing:
#'         \item{n}{number of complete (all energy balance terms available) observations}
#'         \item{intercept}{intercept of the OLS regression}
#'         \item{slope}{slope of the OLS regression}
#'         \item{r_squared}{r^2 of the OLS regression}
#'         \item{EBR}{energy balance ratio}
#' 
#' 
#' @references Wilson K., et al. 2002: Energy balance closure at FLUXNET sites.
#'             Agricultural and Forest Meteorology 113, 223–243.

# include second method: regression (see Wilson 2002) (evtl RMA)
# include case: S not in the dataframe, but not set to NULL in the arguments.

energy.closure <- function(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",
                           missing.G.as.NA=F,missing.S.as.NA=F){
  
  Rn <- check.columns(data,Rn)
  LE <- check.columns(data,LE)
  H  <- check.columns(data,H)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,nrow(data))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,nrow(data))
  }


  comp <- complete.cases(Rn,G,S,LE,H)
  n    <- sum(comp)
  
  EBR <- sum(LE[comp] + H[comp]) / sum(Rn[comp] - G[comp] - S[comp])
  
  emod <- lm((LE + H) ~ (Rn - G - S))
  intercept <- summary(emod)$coef[1,1]  
  slope     <- summary(emod)$coef[2,1] 
  r_squared <- summary(emod)$r.squared
  
  return(c("n"=n,"intercept"=intercept,"slope"=slope,"r^2"=r_squared,"EBR"=EBR))
  
}




## is instanteneous even the right approach here?

#' instanteneous energy balance closure and flux adjustment
#' 
#' @description Characterization and correction of instantaneous energy balance
#'              non-closure. Fluxes are adjusted based on two common methods: 
#'              Bowen ratio and ...
#' 
#' @param data data.frame or matrix containing all required columns
#' @param Rn   net radiation (W m-2)
#' @param G    ground heat flux (W m-2)
#' @param S    sum of all storage fluxes (W m-2)
#' @param LE   latent heat flux (W m-2)
#' @param H    sensible heat flux (W m-2)
#' 
#' 

energy.closure <- function(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",
                           missing.G.as.NA=F,missing.S.as.NA=F){
  
  Rn <- check.columns(data,Rn)
  LE <- check.columns(data,LE)
  H  <- check.columns(data,H)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,nrow(data))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,nrow(data))
  }
  

  comp <- complete.cases(Rn,G,S,LE,H)
  
  energy.gap <- (Rn[comp] - G[comp] - S[comp]) - (LE[comp] + H[comp])
  
  bowen.ratio <- H/LE
  
  H.adjust  <- energy.gap * bowen.ratio

  
}



## General literature
## Gb: Hong_2012

