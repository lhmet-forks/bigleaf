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

### Issues/check:
# MOL independent of temperature?
# more efficient way to provide data.frames? rather than specifiyng each column separately?
# check if Martin 1989, Daudet 1999 etc. really need Gb only or Ga?
# read + understand the differences of the omega formulations, starting with appendix of Jarvis & McNaughton 1986!

## TODO:
# Esat <- provide only one option? Else this option has to be inlcuded in all functions where Esat is used!!
#         not really, if it was  to be changed, do it in the function directly, not where it is a sub-function!
#        options that are not often used may not be indicated: rather changed in the subfunction directly --> ensures consistency! 
# check e to q conversions again!
# References for EVERY function!!!
# check for more functions in Jarvis,McNaughton and the Monteith textbook!!!
# doublecheck Choudhury!
# write examples!
# think about combining Gb functions into one function!!!!!
# virtual temperature. replace q with VPD!
# include option that variables can also be a vector? but not necessarily...maybe in "check.columns"

### still missing:
# ETpot - add other formulations: see Donahue 2010!
# Ga in mol and Ga for CO2!!
# inst. energy balance closure + correction
# Uncertainty!!
# decoupling: Martin_1989 hypostomatous
# dew point temperature
# potential temperature



# Notes:
# other estimates for d and z0m: Shaw and Pereira 1982 (see also Yang_2003)
# v = kinematic viscosity of the air (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html


## artificial test data
# wind     <- rnorm(20,3,0.8)
# ustar    <- rnorm(20,0.5,0.2)
# Tair     <- rnorm(20,25,2)
# pressure <- rep(101.325,20)
# Rn       <- rnorm(20,400,40)
# H        <- 0.3*Rn
# LE       <- 0.4*Rn
# G        <- 0.1*Rn
# Ga       <- rnorm(20,0.08,0.01)
# Gs       <- rnorm(20,0.3,0.02)*(1/41)
# VPD      <- rnorm(20,2,0.3)
# precip   <- c(rep(0,18),0.1,0.001)
# PPFD     <- seq(180,300,length.out=20)
# GPP      <- rnorm(20,20,3)
# NEE      <- -GPP
# day      <- c(rep(1,10),rep(2,10))
# 
# 
# 
# testdf <- data.frame(wind,ustar,Tair,pressure,Rn,
#                      H,LE,G,Ga,Gs,VPD,precip,PPFD,GPP,NEE,day)


## real data
# path.pkg <- "C:/Profiles/jknauer/Documents/RPackage_BigLeaf/"
# load(paste0(path.pkg,"DK-Sor_2010_processed.rda"))
# testdk <- DK_Sor_2010_processed
# colnames(testdk)[10] <- "VPD"
# testdk[,"VPD"] <- testdk[,"VPD"]/10
# save(testdk,file="C:/Profiles/jknauer/Desktop/testdk.RData")


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
    Rbwc       = 1.4              # Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
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
#' @param Ca            atmospheric CO2 concentration (umol mol-1)
#' @param day           day of year
#' @param month         month
#' @param year          year
#' @param quality_ext   
#' @param good_quality 
#' @param excl_vars_qc  
#' @param tprecip       precipitation threshold used to demark a precipitation event (mm)
#' @param precip_hours  number of hours removed following a precipitation event (h)
#' @param trad          radiation threshold (umol m-2 s-1)
#' @param ttemp         temperature threshold (deg C)
#' @param tustar        friction velocity threshold (m s-1)
#' @param tGPP          GPP threshold (???)
#' @param trH           relative humidity threshold (-)
#' @param NA.as.precip  if TRUE missing precipitation data are treated as precipitation events
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

filter.data <- function(data,precip="precip",PPFD="PPFD",Tair="Tair",ustar="ustar",
                        VPD="VPD",GPP="GPP_nt",Ca="Ca",day="day",month="month",year="year",
                        quality_ext="_qc",good_quality=c(0,1),excl_vars_qc=c("PPFD","Ca"),
                        tprecip=0.01,precip_hours=24,trad=200,ttemp=5,tustar=0.2,tGPP=0.5,
                        trH=0.95,NA.as.precip=F){
  
  vars    <- c(precip,PPFD,Tair,ustar,VPD,GPP,Ca)
  vars_qc <- setdiff(vars,excl_vars_qc)
  
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
  
  if (Ca %in% vars_qc){
    Ca   <- check.columns(data,Ca) 
  }
  
  ## data quality
  cat("Quality control:",fill=TRUE)
  for(var in vars_qc){
    if (paste0(var,quality_ext) %in% colnames(data)){                             # does column exist?
      assign(paste0(var,quality_ext),check.columns(data,paste0(var,quality_ext))) # create variable "*_qc", make quality check before
      data[get(paste0(var,quality_ext)) > max(good_quality),var] <- NA            # exclude bad quality data
      
      qc_invalid      <- sum(get(paste0(var,quality_ext)) > max(good_quality)) # count & report
      qc_invalid_perc <- round((qc_invalid/nrow(data))*100,2)
      
      cat(var,": ",qc_invalid," data points (",qc_invalid_perc,"%) excluded",fill=T,sep="")
    } 
  }
  cat("-------------------------------------",fill=T)
  
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
  if (NA.as.precip){
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
  
  addition_precip <- length(setdiff(precip_invalid,unique(growseas_invalid)))
  addition_PPFD   <- length(setdiff(PPFD_invalid,unique(c(growseas_invalid,precip_invalid))))
  addition_Tair   <- length(setdiff(Tair_invalid,unique(c(growseas_invalid,precip_invalid,PPFD_invalid))))
  addition_ustar  <- length(setdiff(ustar_invalid,unique(c(growseas_invalid,precip_invalid,PPFD_invalid,Tair_invalid))))
  addition_rH     <- length(setdiff(rH_invalid,unique(c(growseas_invalid,precip_invalid,PPFD_invalid,Tair_invalid,ustar_invalid))))
  
  addition_precip_perc <- round(addition_precip/nrow(data)*100,2)
  addition_PPFD_perc <- round(addition_PPFD/nrow(data)*100,2)
  addition_Tair_perc <- round(addition_Tair/nrow(data)*100,2)
  addition_ustar_perc <- round(addition_ustar/nrow(data)*100,2)
  addition_rH_perc <- round(addition_rH/nrow(data)*100,2)
  
  cat("Data filters:",fill=T)
  cat(length(growseas_invalid)," data points (",growseas_perc,"%) excluded by growing season filter",fill=T,sep="")
  cat(addition_precip," additional data points (",addition_precip_perc,"%) excluded by precipitation filter (",length(precip_invalid)," data points = ",precip_perc,"% in total)",fill=T,sep="")
  cat(addition_PPFD," additional data points (",addition_PPFD_perc,"%) excluded by radiation filter (",length(PPFD_invalid)," data points = ",PPFD_perc,"% in total)",fill=T,sep="")
  cat(addition_Tair," additional data points (",addition_Tair_perc,"%) excluded by air temperature filter (",length(Tair_invalid)," data points = ",Tair_perc,"% in total)",fill=T,sep="")
  cat(addition_ustar," additional data points (",addition_ustar_perc,"%) excluded by friction velocity filter (",length(ustar_invalid)," data points = ",ustar_perc,"% in total)",fill=T,sep="")
  cat(addition_rH," additional data points (",addition_rH_perc,"%) excluded by relative humidity filter (",length(rH_invalid)," data points = ",rH_perc,"% in total)",fill=T,sep="")
  
  invalid        <- unique(c(growseas_invalid,precip_invalid,PPFD_invalid,Tair_invalid,ustar_invalid,rH_invalid))
  valid[invalid] <- 0
  
  excl_perc <- round((length(invalid)/nrow(data))*100,2)
  
  cat(length(invalid)," data points (",excl_perc,"%) excluded in total",fill=T,sep="")
  cat(nrow(data) - length(invalid)," valid data points (",100-excl_perc,"%) remaining.",fill=T,sep="")
  
  # 5) set all other columns to NA
  data_filtered <- data.frame(data,valid)
  data_filtered[valid==0,!colnames(data) %in% c("day","month","year")] <- NA
  
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
#' @details none 



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
#' @param data      data.frame or matrix containing all required variables
#' @param GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' @param NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Tair      Air temperature (deg C)
#' @param constants 
#'
#' @details the following metrics are calculated:
#' 
#'          Water-use efficiency (WUE):
#'          
#'          \deqn{WUE = GPP / ET}
#'          
#'          Water-use efficiency based on NEE (WUE_NEE):
#'          
#'          \deqn{WUE_NEE = NEE / ET}
#'          
#'          Inherent water-use efficiency (IWUE; Beer et al. 2009):
#'          
#'          \deqn{IWUE = (GPP * VPD) / ET}
#'          
#'          Underlying water-use efficiency (uWUE; Zhou et al. 2014):
#'          
#'          \deqn{uWUE= (GPP * sqrt(VPD)) / ET}
#' 
#' @return a named vector with the following elements:
#'         \item{WUE}{Water-use efficiency ()}
#'         \item{WUE_NEE}{Water-use efficiency based on NEE ()}
#'         \item{IWUE}{Inherent water-use efficiency ()}
#'         \item{uWUE}{Underlying water-use efficiency ()}
#' 
#' @note Units for VPD can also be hPa. Units change accordingly.
#' 
#' @references Beer, C., et al., 2009: Temporal and among-site variability of inherent
#'             water use efficiency at the ecosystem level. Global Biogeochemical Cycles 23, GB2018.
#'             
#'             Zhou, S., et al., 2014: The effect of vapor pressure deficit on water
#'             use efficiency at the subdaily time scale. Geophysical Research Letters 41.
#'                                      
#' @export
WUE.metrics <- function(data,GPP="GPP_nt",NEE="NEE",LE="LE",VPD="VPD",Tair="Tair",
                        constants=bigleaf.constants()){
  
  GPP  <- check.columns(data,GPP)
  NEE  <- check.columns(data,NEE)
  LE   <- check.columns(data,LE)
  VPD  <- check.columns(data,VPD)
  Tair <- check.columns(data,Tair)
  
  
  ET  <- LE.to.ET(LE,Tair)                   # kg H2O s-1
  GPP <- (GPP/1e06 * constants$Cmol)*1000  # gC m-2 s-1
  NEE <- (NEE/1e06 * constants$Cmol)*1000  # gC m-2 s-1
  
  
  WUE     <- median(GPP/ET,na.rm=T)
  WUE_NEE <- median(abs(NEE)/ET,na.rm=T)
  IWUE    <- median((GPP*VPD)/ET,na.rm=T)
  uWUE    <- median((GPP*sqrt(VPD))/ET,na.rm=T)
  
  return(c(WUE=WUE,WUE_NEE=WUE_NEE,IWUE=IWUE,uWUE=uWUE))
}


# by(data=testdk2,testdk2[,"month"],WUE.metrics,simplify=T) -> test
# do.call(cbind,test) -> test2  # convert to data.frame


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


################################################
### Boundary layer conductance formulations #### -------------------------------------------------------------------------
################################################

#' Boundary layer conductance according to Thom 1972
#' 
#' @description An empirical formulation for the canopy boundary layer conductance
#'              based on a simple ustar dependency.
#' 
#' @param ustar     friction velocity (m s-1)
#' @param constants k - von-Karman constant (-) \cr
#'                  Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#' 
#' @return a data.frame with the following columns:
#'  \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'  \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
#'  \item{Gb}{Boundary layer conductance (m s-1)}
#'  \item{kB}{kB-1 parameter (-)}
#'  
#'  
#' @details
#'  
#'  The empirical equation for Rb to water suggested by Thom 1972 is:
#'  
#'  \deqn{Rb = 6.2ustar^-0.67}
#'  
#'  Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'  \deqn{Rb_CO2 = 1.4 * Rb}
#'  
#'  The factor 1.4 arises due the lower molecular diffusivity of CO2 compared to water.
#'  It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'  across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'  mixing (Nobel 2005).
#' 
#' @references Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
#'             Quarterly Journal of the Royal Meteorological Society 98, 124-134.
#'             
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#' 
#' @seealso \code{\link{Gb.Su}}, \code{\link{Gb.Choudhury}}
#' 
#' 
#' @export
Gb.Thom <- function(ustar,constants=bigleaf.constants()){
  Rb     <- 6.2*ustar^-0.667
  Gb     <- 1/Rb
  kB     <- Rb*constants$k*ustar
  Rb_CO2 <- constants$Rbwc * Rb
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}


#' Roughness Reynolds Number
#' 
#' @description 
#' 
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param hs        Roughness length of the soil (m)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K)
#'                  
#' @details The roughness Reynolds Number is calculated as in Massman 1999a:
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





#' Boundary layer conductance according to Su et al. 2001
#' 
#' A physically based formulation for the canopy boundary layer conductance. 
#'
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      air temperature (degC)
#' @param presssure atmospheric pressure (Pa)
#' @param ustar     friction velocity (m s-1)
#' @param wind      wind speed (m s-1)
#' @param H         sensible heat flux (W m-2)
#' @param zh        canopy height (m)
#' @param zr        reference height (m)
#' @param d         zero-plane displacement height (m)
#' @param z0m       roughness length for momentum (m)
#' @param Dl        leaf dimension (m)
#' @param N         number of leaf sides participating in heat exchange (1 or 2)
#' @param fc        fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI       one-sided leaf area index (-)
#' @param Cd        foliage drag coefficient (-)
#' @param hs        roughness height of the soil (m)
#' @param stab_correction should stability correction be applied? Defaults to TRUE
#' @param stab_formulation stability correction function used (If stab_correction is TRUE).
#'                         Either "Dyer_1970" or "Businger_1971".
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K) \cr
#'                  Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#' 
#' @return a data.frame with the following columns:
#'     \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'     \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
#'     \item{Gb}{Boundary layer conductance (m s-1)}
#'     \item{kB}{kB-1 parameter (-)}
#'     
#' @details The formulation is based on the kB-1 model developed by Massman 1999. 
#'          Su et al. 2001 derived the following approximation:
#'           
#'          \deqn{kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2}
#'          
#'          If fc (fractional vegetation cover) is missing, it is estimated from LAI:
#' 
#'          \deqn{fc = 1 - exp(-LAI/2)}
#'          
#'          The wind speed at the top of the canopy is calculated using function
#'          \code{\link{wind.profile}}.
#'          
#'          Ct is the heat transfer coefficient of the leaf (Massman 1999):
#'          
#'          \deqn{Ct = Pr^-2/3 Reh^-1/2 N}
#'          
#'          where Pr is the Prandtl number (set to 0.71), and Reh is the Reynolds number for leaves:
#'          
#'          \deqn{Dl wind(zh) / v}
#'           
#'          kBs-1, the kB-1 value for bare soil surface, is calculated according 
#'          to Su et al. 2001:
#'          
#'          \deqn{kBs^-1 = 2.46(Re)^0.25 - ln(7.4)}
#'          
#'          Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'          \deqn{Rb_CO2 = 1.4 * Rb}
#'  
#'          The factor 1.4 arises due the lower molecular diffusivity of CO2 compared to water.
#'          It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'          across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'          mixing (Nobel 2005).
#' 
#' @references Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
#'             two models for estimation of the roughness height for heat transfer between
#'             the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
#' 
#'             Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#'            
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#' 
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Choudhury}}
#' 
#' @export
Gb.Su <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",wind="wind",
                  H="H",zh,zr,d,z0m,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,stab_correction=TRUE,
                  stab_formulation=c("Dyer_1970","Businger_1971"),
                  constants=bigleaf.constants()){
        
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  wind     <- check.columns(data,wind)
  H        <- check.columns(data,H)
  
  if (is.null(fc)) {
    fc  <- (1-exp(-LAI/2)) 
  } 
  
  wind_zh <- wind.profile(zh,data,zr=zr,d=d,z0m=z0m,stab_correction=stab_correction,
                          stab_formulation=stab_formulation)[,1]
  
  v   <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"v"]
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"Re"]
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind_zh / v
  Ct  <- 1*0.71^-0.6667*Reh^-0.5*N   # 0.71 = Prandtl number                       

  kB     <- (constants$k*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb     <- kB/(constants$k*ustar) 
  Gb     <- 1/Rb
  Rb_CO2 <- constants$Rbwc * Rb
  
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}

## typical values for hs: Shuttleworth and Wallace 1985 from van Bavel and Hillel 1976 -> 0.01m

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
#' @param wind             Wind speed at sensor height (m s-1)
#' @param ustar            Friction velocity (m s-1)
#' @param leafwidth        Leaf width (m)
#' @param LAI              One-sided leaf area index
#' @param zh               Canopy height (m)
#' @param zr               Instrument (reference) height (m)
#' @param stab_correction  Should stability correction be applied? Defaults to TRUE
#' @param stab_formulation Stability correction function used (If stab_correction is TRUE).
#'                         Either "Dyer_1970" or "Businger_1971".
#' @param constants        k - von-Karman constant (-) \cr
#'                         Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#' 
#' @return a data frame with the following columns:
#'     \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'     \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
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
#'          Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'          \deqn{Rb_CO2 = 1.4 * Rb}
#'  
#'          The factor 1.4 arises due the lower molecular diffusivity of CO2 compared to water.
#'          It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'          across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'          mixing (Nobel 2005).
#'          
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             McNaughton, K. G., Van den Hurk, B.J.J.M., 1995: A 'Lagrangian' revision of
#'             the resistors in the two-layer model for calculating the energy budget of a
#'             plant canopy. Boundary-Layer Meteorology 74, 261-288.
#'             
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#'             
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Su}}  
#'            
#' @export                                                                                                                                                                                                                                                                                    
Gb.Choudhury <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                         leafwidth,LAI,zh,zr,d,z0m,stab_correction=TRUE,
                         stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  ustar    <- check.columns(data,ustar)
  wind     <- check.columns(data,wind)
  
  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)
  wind_zh <- wind.profile(zh,data,zr=zr,d=d,z0m=z0m,stab_correction=stab_correction,
                          stab_formulation=stab_formulation)[,1]
  Gb      <- LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  Rb      <- 1/Gb
  kB      <- Rb*constants$k*ustar
  Rb_CO2  <- constants$Rbwc * Rb
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}


#' Wind speed at given levels in the surface layer
#' 
#' @description wind speed at a given height above the canopy estimated from single-level
#'              measurements of wind speed at some distance above the canopy.
#'              
#' @details The underlying assumption is the existence of a logarithmic wind profile
#'          above the height d + z0m (the height at which wind speed mathematically reaches zero
#'          according to the Monin-Obhukov similarity theory).
#'          In this case, the wind speed at a given height z is given by:
#'          
#'          \deqn{u(z) = (ustar/k) * (ln((z - d) / z0m) - \psi{m}}
#'          
#'          The roughness parameters zero-plane displacement height (d) and roughness length (z0m)
#'          can be approximated from \code{\link{roughness.parameters}}.
#'          
#' @param heights   vector with heights for which wind speed is to be 
#'                  calculated.
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)                                                                                  
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zr        Instrument (reference) height (m)
#' @param d         Zero-plane displacement height (m)
#' @param z0m       Roughness length for momentum (m)
#' @param stab_correction Should stability correction be applied? Defaults to TRUE
#' @param stab_formulation Stability correction function used (If stab_correction is TRUE).
#'                         Either "Dyer_1970" or "Businger_1971".
#' @param constants k - von-Karman constant (-) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  g - gravitational acceleration (m s-2) \cr
#'                                   
#' @note Note that this equation is only valid for z >= d + z0m. All values in \code{heights}
#'       smaller than d + z0m will return 0.                                 
#'                                  
#' @return a data.frame with rows representing time and columns representing heights.     
#'                                            
#' @references
#' 
#' @export                                                                                                                          
wind.profile <- function(heights,data,Tair="Tair",pressure="pressure",ustar="ustar",
                         H="H",zr,d,z0m,stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  if( any(heights < d + z0m) ) warning("function is only valid for heights above d + z0m! Wind speed for heights below d + z0m will return 0!") 
  
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  
  wind_heights <- data.frame(matrix(NA,ncol=length(heights),nrow=length(Tair)))
  colnames(wind_heights) <- paste0(heights,"m")
  for (z in heights){
    i <- which(heights == z)
    
    if (stab_correction){
      psi_m <- stability.correction(data=data,zr=z,d=d,
                                    formulation=stab_formulation,constants=constants)[,"psi_m"]
      
      wind_heights[,i] <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m) - psi_m))
    
    } else {
      
      wind_heights[,i] <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m)))
      
    }
    
  }
  
  return(wind_heights)
}



#' Roughness parameters
#' 
#' @description A simple approximation of the two roughness parameters displacement height (d)
#'              and roughness length for momentum (z0m).
#'              
#' @param method    Method used, either "canopy_height", or "wind_profile" \cr
#'                  NOTE: if method is "canopy_height", only the following three arguments
#'                  are used.               
#' @param frac_d    Fraction of displacement height on canopy height (-)
#' @param frac_z0m  Fraction of roughness length on canopy height (-)
#'                  The following arguments are only needed if method = "wind_profile"!
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param wind      Wind speed at height zr (m s-1)
#' @param ustar     Friction velocity (m s-1)
#' @param H              Sensible heat flux (W m-2)
#' @param zr             Instrument (reference) height (m)
#' @param d              Zero-plane displacement height (m); optional
#' @param z0m            Roughness length for momentum (m); optional
#' @param stab_roughness   Should stability correction be considered? Default is TRUE
#' @param stab_formulation Stability correction function used (If stab_correction is TRUE).
#'                         Either "Dyer_1970" or "Businger_1971".
#' @param constants k - von-Karman constant (-) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  g - gravitational acceleration (m s-2) \cr
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
#'          \deqn{z0m = median((zr - d) * exp(-k*wind / ustar - psi_m)}
#'                  
#'          By default, d in this equation is fixed to 0.7zh, but can be set to any
#'          other value.        
#' 
#' @return a data.frame with the following columns:
#'         \item{d}{Zero-plane displacement height (m)}
#'         \item{z0m}{Roughness length for momentum (m)}
#'
#'    
#' @return                                  
roughness.parameters <- function(method=c("wind_profile","canopy_height"),zh,
                                 frac_d=0.7,frac_z0m=0.1,data,Tair="Tair",pressure="pressure",
                                 wind="wind",ustar="ustar",H="H",zr,d=NULL,z0m=NULL,
                                 stab_roughness=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf.constants()){
  
  method           <- match.arg(method)
  stab_formulation <- match.arg(stab_formulation)
  
  if (method == "canopy_height"){
  
    d      <- frac_d*zh
    z0m    <- frac_z0m*zh
    z0m_se <- NA
  
  } else if (method == "wind_profile"){
    
    Tair     <- check.columns(data,Tair)
    pressure <- check.columns(data,pressure)
    wind     <- check.columns(data,wind)
    ustar    <- check.columns(data,ustar)
    H        <- check.columns(data,H)
    
    
    if (is.null(d) & is.null(z0m)){
      
      d <- frac_d*zh
    
    }
    
    if (stab_roughness){
      
      psi_m <- stability.correction(data=data,zr=zr,d=d,
                                    formulation=stab_formulation,constants=constants)[,"psi_m"]

      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar - psi_m)
      
    } else {
      
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar)
      
    }
    
    z0m_all[z0m_all > 1000] <- NA
    
    z0m    <- median(z0m_all,na.rm=T)
    z0m_se <- 1.253 * (sd(z0m_all,na.rm=T) / sqrt(sum(!is.na(z0m_all))))

  }

  return(data.frame(d,z0m,z0m_se))
}

# psi_h <- stability.correction(data,Tair=testdk2[,"Tair"],pressure=testdk2[,"pressure"],ustar=testdk2[,"ustar"],
#                               H=testdk2[,"H"],zr=43,d=0.7*25,)[,"psi_m"]
# 
# 
# ram1 <- wind/ustar^2
# ram2 <- log((43-17.5)/1.65)/(k*ustar)
# 
# ram3 <- log((43-17.5)/1.07 - psi_h)/(k*ustar)
# ram4 <- log((43-17.5)/2.5 - psi_h)/(k*ustar)

############################
### air pressure, density ## --------------------------------------------------------------------------
############################

#' Air density
#' 
#' @description Air density of moist air from air temperature and pressure.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param constants Kelvin - conversion degC to Kelvin \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1)
#' 
#' @details Air density (\eqn{\rho}) is calculated as:
#' \deqn{\rho = pressure / Rd * Tair}
#' 
#' @return rho air density (kg m-3)
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
#' @param elev      Elevation a.s.l. (m)
#' @param Tair      Air temperature (deg C)
#' @param q         Specific humidity (kg kg-1); optional
#' @param constants Kelvin- conversion degC to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  g -  gravitational acceleration (m s-2) \cr
#' 
#' 
#' @details 
#' 
#' @note The hypsometric equation gives an estimate of the standard pressure
#'       at a given altitude. 
#'       If specific humidity q is provided, humidity correction
#'       is applied and the virtual temperature instead of air temperature is used.
#'
#' @return pressure air pressure (kPa)
#'                            
#' @references Stull B., 1988: An Introduction to Boundary Layer Meteorology.
#'             Kluwer Academic Publishers, Dordrecht, Netherlands.
#' 
#' @examples
#' # mean pressure at 500m altitude at 25 deg C and specific humidity of 0.01 kg kg-1
#' pressure.from.elevation(500,Tair=25,q=0.01)
#' 
#' @export                           
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
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
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
#' @note Note that L gets very small for very low ustar values with implications
#'       for subsequent functions using L as input. It is recommended to filter
#'       data and exclude low ustar values (ustar < ~0.2) beforehand. 
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
#' @param Tair      Air temperature (degC)
#' @param pressure  Air pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zr        Instrument (reference) height (m)
#' @param disp      Zero-plane displacement height (m)
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


#' stability correction functions for heat and momentum
#' 
#' gives the integrated form of the universal functions.
#' 
#' @description dimensionless stability functions needed to correct deviations
#'              from the exponential wind profile under non-neutral conditions.
#'              
#' @param Tair         Air temperature (deg C)
#' @param pressure     Air pressure (kPa)
#' @param ustar        Friction velocity (m s-1)
#' @param H            Sensible heat flux (W m-2)
#' @param zr           instrument (reference) height (m)
#' @param d            Zero-plane displacement height (m)
#' @param formulation  Formulation for the stability function. Either "Dyer_1970", or "Businger_1971"
#' @param constants    Kelvin - conversion degree Celsius to Kelvin \cr
#'                     Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                     cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                     k - von Karman constant (-) \cr
#'                     g - gravitational acceleration (m s-2)  
#'
#' @details The functions depend on the value of the stability parameter \eqn{\zeta},
#'          which can be calculated from the function \code{\link{stability.parameter}}.
#'          The integration of the universal functions is:
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
#' @return a data.frame with the following columns:
#'          \item{psi_h}{the value of the stability function for heat and water vapor (-)}
#'          \item{psi_m}{the value of the stability function for momentum (-)}
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
#'             Paulson, C.A., 1970: The mathematical representation of wind speed
#'             and temperature profiles in the unstable atmospheric surface layer.
#'             Journal of Applied Meteorology 9, 857-861.
#' 
#'             Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#'   
stability.correction <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",H="H",zr,d,
                                 formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf.constants()){
  
  formulation  <- match.arg(formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  
  zeta <- stability.parameter(Tair,pressure,ustar,H,zr,d,constants)
  
  psi_h = psi_m <- numeric()
  
  # universal functions
  if (formulation == "Businger_1971"){
    x_h <- -7.8
    x_m <- -6
    y_h <- 0.95 * ( 1 - 11.6 * zeta)^0.5
    y_m <- (1 - 19.3*zeta)^0.25
  } else if (formulation == "Dyer_1970"){
    x_h = x_m <- -5
    y_h       <- (1 - 16 * zeta)^0.5
    y_m       <- (1 - 16 * zeta)^0.25
  }

  # integration of universal functions (after Paulson_1970 and Foken 2008)
  # stable
  stable <- zeta >= 0 | is.na(zeta)
  psi_h[stable] <- x_h * zeta[stable]
  psi_m[stable] <- x_m * zeta[stable]
  # unstable
  unstable <- zeta < 0 | is.na(zeta)
  psi_h[unstable] <- 2 * log( (1 + y_h[unstable] ) / 2)  
  psi_m[unstable] <- 2 * log( (1 + y_m[unstable] ) / 2) + 
                     log( ( 1 + y_m[unstable]^2 ) / 2) 
                     -2 * atan(y_m[unstable]) + pi/2
  
  return(data.frame(psi_h,psi_m))
} 

# zeta <- seq(-5,-0.5,0.5)
# x1<- 2*log(( 1 + (1 - 16 * zeta)^0.5) / 2)
# x2<- 2*log(( 1 + 1/(1 - 16 * zeta)^-0.5) / 2)
# 
# y1 <- 2 * log( 1 +   (0.95*(1 - 11.6*zeta)^0.5) / 2)
# y2 <- 2 * log( 1 + 1/(0.95*(1 - 11.6*zeta)^-0.5) / 2)
# 
# xx1 <- 0.74*(1 - 9*zeta)^-0.5
# xx2 <- 0.74 + 4.7*zeta
### Dyer_1970; Dyer_1974 (similar:Paulson 1970)
### Businger_1971: from Foken_2008; p.65 (suggested by Businger_1971), in the form of Högstrom 1988



#' Aerodynamic conductance
#' 
#' @description Aerodynamic conductance, including options for the boundary layer conductance
#'              formulation and stability correction functions.
#' 
#' @param data              Data.frame or matrix containing all required variables
#' @param Tair              Air temperature (deg C)
#' @param pressure          Air pressure (kPa)
#' @param wind              Wind speed (m s-1)
#' @param ustar             Friction velocity (m s-1)
#' @param H                 Sensible heat flux (W m-2)
#' @param zr                Instrument (reference) height (m)
#' @param zh                Canopy height (m)
#' @param d                 Zero-plane displacement height (m)
#' @param z0m               Roughness length for momentum (m)
#' @param Dl                Characteristic leaf dimension (m); only used if \code{Rb_model} is "Choudhury_1988" or "Su_2001".
#' @param N                 Number of leaf sides participating in heat exchange (1 or 2); only used if \code{Rb_model} is "Su_2001".
#' @param fc                Fractional vegetation cover (-); only used if \code{Rb_model} is "Su_2001". See Details.
#' @param LAI               One-sided leaf area index (m2 m-2); only used if \code{Rb_model} is "Choudhury_1988" or "Su_2001".
#' @param Cd                Foliage drag coefficient (-); only used if \code{Rb_model} is "Su_2001". 
#' @param hs                Roughness length of bare soil (m); only used if \code{Rb_model} is "Su_2001".
#' @param stab_correction   Should stability correction be applied? Default is TRUE
#' @param stab_formulation  Stability correction function. Either "Dyer_1970" (default) or
#'                          "Businger_1971". Ignored if \code{stab_correction} is FALSE.
#' @param Rb_model          Boundary layer resistance formulation. One of c("Thom_1972","Choudhury_1988","Su_2001).
#' @param kB                kB-1 value; only used if \code{Rb_model = "constant_kB"}
#' @param constants         k - von Karman constant (-) \cr
#'                          Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-) \cr
#'                          cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                          Kelvin - conversion degree Celsius to Kelvin \cr
#'                          g - gravitational acceleration (m s-2) \cr
#'                          pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                          Tair0 - reference air temperature (K) \cr
#'
#' @details
#' 
#' Aerodynamic conductance (Ga) is calculated as:
#' 
#'  Ga = 1 / (Ra_m + Rb)
#'  
#'  where Ra_m is the aerodynamic resistance for momentum and Rb the canopy boundary
#'  layer resistance.
#'  
#'  The argument \code{stab_formulation} determines the stability correction function used 
#'  to account for the effect of atmospheric stability on Ra_m (Ra_m is lower for unstable
#'  and higher for stable stratification). Stratification is based on a stability parameter zeta (z-d/L),
#'  where z = reference height, d the zero-plane displacement height, and L the Monin-Obukhov length, 
#'  calculated with \code{\link{MoninObukhov.length}}
#'  Stability correction functions calculate the correction based on the formulations by Businger 1971 
#'  Dyer 1970, based on zeta.
#'  If \code{stab_correction} is set to FALSE, the effects of atmospheric stability on Ga are neglected,
#'  and the Monin-Obukhov length, the stability parameter zeta, as well as the stability correction functions
#'  are not calculated.
#' 
#'  The argument \code{Rb_model} determines the canopy boundary layer resistance ("excess resistance") model that is used. "Thom_1972" is an 
#'  empirical formulation based on the friction velocity (ustar) (Thom 1972):
#'  
#'    \deqn{Rb = 0.667u*^0.667}
#'    
#'  The model by Choudhury 1988 ("Choudhury_1988"), calculates Rb
#'  based on leaf width, LAI and ustar (Note that in the original formulation leaf width
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
#'  Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'  \deqn{Rb_CO2 = 1.4 * Rb}
#'  
#'  The factor 1.4 arises due the lower molecular diffusivity of CO2 compared to water.
#'  It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'  across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'  mixing (Nobel 2005).
#'  
#'  If the roughness parameters z0m and d are not available, they can be estimated using
#'  \link{\code{roughness.parameters}}.
#'  
#' @return a dataframe with the following columns:
#'         \item{Ga_m}{Aerodynamic conductance for momentum (m s-1)}
#'         \item{Ga_h}{Aerodynamic conductance for heat and water vapor (m s-1)}
#'         \item{Ga_CO2}{Aerodynamic conductance for CO2 (m s-1)}
#'         \item{Ra_m}{Aerodynamic resistance for momentum (s m-1)}
#'         \item{Ra_h}{Aerodynamic resistance for heat and water vapor (s m-1)}
#'         \item{Ra_CO2}{Aerodynamic resistance for CO2 (s m-1)}
#'         \item{Gb}{Canopy boundary layer conductance for heat and water vapor (m s-1)}
#'         \item{Rb}{Canopy boundary layer resistance for heat and water vapor (s m-1)}
#'         \item{Rb_CO2}{Canopy boundary layer resistance for CO2 (s m-1)}
#'         \item{kB}{kB-1 parameter (-)}
#'         \item{zeta}{Stability parameter 'zeta' (-)}
#'         \item{psi_h}{Integrated stability correction function (-)}
#'        
#'         
#' @references 
#'  
#' @export
#'    


# TO DO:
# - find out what "." means in front of a function!)
# - think about how to deal with uncertainty!
# - think about testing for missing arguments in the subfunctions separately, but most likely doesn't make sense
# - better delineate the optionality of the arguments!
aerodynamic.conductance <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                                    zr,zh,d,z0m,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,
                                    stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                    Rb_model=c("Thom_1972","Choudhury_1988","Su_2001","constant_kB"),
                                    kB=NULL,constants=bigleaf.constants()){
  
  Rb_model         <- match.arg(Rb_model)
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  wind     <- check.columns(data,wind)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  

  ## calculate boundary layer conductance
  if (Rb_model != "constant_kB"){
    
    if (Rb_model == "Thom_1972"){
      
      Gb_mod <- Gb.Thom(ustar=ustar,constants=constants)
      
    } else if (Rb_model == "Choudhury_1988"){
      
      Gb_mod <- Gb.Choudhury(data,leafwidth=Dl,LAI=LAI,zh=zh,zr=zr,d=d,z0m=z0m,
                             stab_correction=stab_correction,
                             stab_formulation=stab_formulation,constants=constants)
      
    } else if (Rb_model == "Su_2001"){
      
      Gb_mod <- Gb.Su(data=data,zh=zh,zr=zr,d=d,z0m=z0m,Dl=Dl,N=N,fc=fc,LAI=LAI,
                      Cd=Cd,hs=hs,stab_correction=stab_correction,
                      stab_formulation=stab_formulation,constants=constants)  
      
    }
 
      kB     <- Gb_mod[,"kB"]
      Rb     <- Gb_mod[,"Rb"]
      Rb_CO2 <- Gb_mod[,"Rb_CO2"]
      Gb     <- Gb_mod[,"Gb"]
  
  } else {
    
    if(is.null(kB)){
      stop("value of kB has to be specified if Rb_model is set to 'constant_kB'!")
    } else {
      Rb     <- kB / (constants$k * ustar)
      Rb_CO2 <- constants$Rbwc * Rb
      Gb     <- 1 / Rb
    }
  
  }
  
  Ra_m <- wind/ustar^2
  
  
  if (stab_correction){
    
    zeta  <-  stability.parameter(Tair,pressure,ustar,H,zr,d,constants)
    
    
    if (stab_formulation %in% c("Dyer_1970","Businger_1971")){
      
      psi_h <- stability.correction(data=data,zr=zr,d=d,
                                    formulation=stab_formulation,constants=constants)[,"psi_h"]
      
    }  
    
    Ra_m  <- pmax((log((zr - d)/z0m) - psi_h),1e-10)/(constants$k*ustar)
  
  
  } else {
    
    zeta = psi_h <- rep(NA,length=length(Ra_m))
    
  }
  
  Ga_m   <- 1/Ra_m
  Ra_h   <- Ra_m + Rb
  Ra_CO2 <- Ra_m + Rb_CO2
  Ga_h   <- 1/Ra_h
  Ga_CO2 <- 1/Ra_CO2

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Ga_CO2,Ra_CO2,Rb,Rb_CO2,Gb,kB,zeta,psi_h))
  
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
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
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
#' @details Surface conductance (Gs) is calculated from the inverted Penman-Monteith equation:
#' 
#'  \deqn{Gs = ( LE * Ga * \gamma ) / ( \delta * (Rn-G-S) + \rho * cp * Ga * VPD - LE * ( \delta + \gamma ) )}
#'  
#'  Note that Gs > Gc (canopy conductance) under conditions when a significant fraction of 
#'  ET comes from interception or soil evaporation. 
#'  Available energy A is defined as A = Rn - G - S. If G and/or S are not provided, then A = Rn.
#'  
#'  By default, any missing data in G and S are set to 0. If arguments \code{missing.S.as.NA} or \code{missing.S.as.NA} are 
#'  set to TRUE, Gs will give NA for these timesteps.
#'  
#'  If pressure is not available, it can be approximated by elevation using the 
#'  function \code{\link{pressure.from.elevation}}
#'  
#' 
#' @return a dataframe with the following columns: 
#'  \item{Gs_ms}{Surface conductance in m s-1}
#'  \item{Gs_mol}{Surface conductance in mol m-2 s-1}
#' 
#' @export
surface.conductance <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga",missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                constants=bigleaf.constants()){
   
   Tair     <- check.columns(data,Tair)
   pressure <- check.columns(data,pressure)
   Rn       <- check.columns(data,Rn)
   VPD      <- check.columns(data,VPD)
   LE       <- check.columns(data,LE)
  
   if(!is.null(G)){
     G <- check.columns(data,G)
     if (!missing.G.as.NA){G[is.na(G)] <- 0}
   } else {
     warning("Ground heat flux G is not provided and set to 0.")
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
   Gs_mol <- ms.to.mol(Gs_ms,Tair,pressure) 

   return(data.frame(Gs_ms,Gs_mol))

}



#' Big-leaf surface conditions
#' 
#' @description Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param H         Sensible heat flux (W m-2)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance for heat and water vapor (m s-1)
#' @param calc.Ca   calculate surface CO2 concentration?
#' @param Ca        Atmospheric CO2 concentration (mol mol-1). Required if calc.Ca = TRUE
#' @param NEE       Net ecosystem exchange (umol m-2 s-1). Required if calc.Ca = TRUE
#' @param Ga_CO2    Aerodynamic conductance for CO2 (mol m-2 s-1). Required if calc.Ca = TRUE          
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
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
#'          see \link{\code{aerodynamic.conductance}}). If Ga is replaced by Ga_m (i.e.
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
#'                                       
#' @export 
#' 

# - add Ga and Ga to data.frame??
bigleaf.surface <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",Ca="Ca",
                            VPD="VPD",Ga="Ga",calc.Ca=F,Ga_CO2="Ga_CO2",NEE="NEE",
                            constants=bigleaf.constants()){

  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  LE       <- check.columns(data,LE)
  H        <- check.columns(data,H)
  VPD      <- check.columns(data,VPD)
  Ga       <- check.columns(data,Ga)
  Ga_CO2   <- check.columns(data,Ga_CO2)
  
  if (calc.Ca){
    Ca       <- check.columns(data,Ca)
    NEE      <- check.columns(data,NEE)
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
  Ca_surf <- as.numeric(rep(NA,length(esat)))
  if (calc.Ca){
    Ca_surf <- Ca.surface(Ca,NEE,Ga_CO2)
  }
  
  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
}



#' Canopy surface CO2 concentration
#'
#' @description the CO2 concentration at the canopy surface derived from net ecosystem
#'              CO2 exchange and measured atmospheric CO2 concentration.
#'              
#' @param Ca     atmospheric CO2 concentration (umol mol-1)
#' @param NEE    net ecosystem exchange (umol CO2 m-2 s-1)
#' @param Ga_CO2 aerodynamic conductance for CO2 (mol m-2 s-1)
#' 
#' @note the following sign convention is employed: negative values of NEE denote
#'       net CO2 uptake by the ecosystem.
#' 
#' @return Ca_surf CO2 concentration at the canopy surface (umol mol-1)
#' 
Ca.surface <- function(Ca,NEE,Ga_CO2){
  
  Ca_surf <- Ca + NEE/Ga_CO2
  
  return(Ca_surf)
  
}



#' Decoupling coefficient
#' 
#' @description The canopy-atmosphere decoupling coefficient 'Omega'. 
#' 
#' @param data        Data.frame or matrix containing all required input variables
#' @param Tair        Air temperature (deg C)
#' @param pressure    Air pressure (kPa)
#' @param Ga          Aerodynamic conductance (m s-1)
#' @param Gs          Surface conductance (m s-1)
#' @param approach    Approach used to calculate omega. Either "JarvisMcNaughton_1986" (default)
#'                    or "Martin_1989".
#' @param canopy.type Either "amphistomatous" (default) or "hypostomatous". See Details. 
#' @param LAI         Leaf area index (m2 m-2), only used if approach = "Martin_1989"
#' @param constants   Kelvin - conversion degree Celsius to Kelvin \cr
#'                    cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                    eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                    sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#' 
#' @details The decoupling coefficient Omega ranges from 0 to 1 and quantifies the
#'          linkage of the conditions (foremost humidity and temperature) at the canopy surface
#'          to the ambient air. Values close to 0 indicate well coupled conditions
#'          characterized by high physiological (i.e. stomatal) control on transpiration
#'          and similar conditions at the canopy surface compared to the atmosphere above
#'          the canopy. Values close to 1 indicate the opposite, i.e. decoupled conditions and 
#'          a low stomatal control on transpiration (Jarvis & McNaughton 1986). \cr
#'          The "JarvisMcNaughton_1986" approach (default option) is the original
#'          formulation for the decoupling coefficient, given by (for an amphistomatous 
#'          canopy (i.e. stomata located on both sides of leaves)):
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1}{\epsilon + 1 + \frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 1) / ( \epsilon + 1 + Ga/Gc)}
#'          
#'          where \eqn{\epsilon = \frac{s}{\gamma}}{\epsilon = s/\gamma} is a dimensionless coefficient
#'          with s being the slope of the saturation vapor pressure curve (Pa K-1), and \eqn{\gamma} the 
#'          psychrometric constant (Pa K-1).
#'          
#'          This formulation can be modified for for hypostomatous (stomata located 
#'          on only one side of the leaves) vegetation (e.g. Daudet et al. 1999): 
#'          
#'          \deqn{\Omega = \frac{\epsilon + 2}{\epsilon + 2 + 2\frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 2) / (\epsilon + 2 + 2Ga/Gc)}
#'          
#'          The approach "Martin_1989" by Martin 1989 additionally takes radiative coupling
#'          into account. For amphistomatous vegetation:
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1 + \frac{Gr}{Ga}}{\epsilon + (1 + \frac{Ga}{Gs}) (1 + \frac{Gr}{Ga})}}{%
#'          \Omega = (\epsilon + 1 + Gr/Ga) / (\epsilon + (1 + Ga/Gs) (1 + Gr/Ga))}
#' 
#' @return omega The decoupling coefficient omega (-)
#' 
#' @references Goldberg V., Bernhofer C., 2008: Testing different decoupling coefficients
#'             with measurements and models of contrasting canopies and soil water conditions.
#'             Annales Geophysicae 26, 1977-1992.
#' 
#'             Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49. 
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
                       approach=c("JarvisMcNaughton_1986","Martin_1989"),
                       canopy.type=c("amphistomatous","hypostomatous"),
                       LAI,constants=bigleaf.constants()){

  approach    <- match.arg(approach)
  canopy.type <- match.arg(canopy.type)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  Ga       <- check.columns(data,Ga)
  Gs       <- check.columns(data,Gs)

  delta   <- Esat(Tair)[,"delta"]
  gamma   <- psychrometric.constant(Tair,pressure,constants=constants)
  epsilon <- delta/gamma
  
  if (approach == "JarvisMcNaughton_1986"){
    
    if (canopy.type == "amphistomatous"){
      
      omega <- (epsilon + 1) / (epsilon + 1 + Ga/Gs)
    
    } else if (canopy.type == "hypostomatous"){
    
      omega <- (epsilon + 2) / (epsilon + 2 + 2*Ga/Gs)
      
    }
    
  } else if (approach == "Martin_1989") {
  
    if (is.null(LAI)){
      
      stop("LAI is not provided!")
      
    } else {
      
      Gr    <- Gr.longwave(Tair,LAI,constants=constants)
      
      if (canopy.type == "amphistomatous"){
      
        omega <- (epsilon + 1 + Gr/Ga) / (epsilon + 1 + Ga/Gs + Gr/Gs + Gr/Ga)
      
      } else if (canopy.type == "hypostomatous"){
        
        stop("option not yet implemented...")
      
      }
    }
  }
  
  return(omega)

}


#' longwave radiative transfer conductance of the canopy
#' 
#' @param Tair      Air temperature (deg C)
#' @param LAI       Leaf area index (m2 m-2)
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


### relate more to surface temperature?
#' Dew point temperature
#' 
#' @description calculates the dew point, the temperature to which air must be 
#'              cooled to become saturated (i.e. e = esat(Td))
#'              
#' @param VPD Vapor pressure deficit (kPa)              



#' Potential temperature
#' 
#' 
#' 

### evtl replace q with VPD!
#' Virtual temperature
#' 
#' @description Virtual temperature, defined as the temperature at which dry air would have the same
#'              density as moist air at its actual temperature.
#' 
#' @param Tair      Air temperature (deg C)
#' @param q         Specific humidity (kg kg-1)
#' @param constants Kelvin - conversion degree Celsius to Kelvin
#' 
#' @details the virtual temperature is given by:
#'  
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
   Tv   <- Tair*(1+0.61*q)
   Tv   <- Tv - constants$Kelvin
   return(Tv)
}


#' Saturation vapor pressure (Esat) and slope of the Esat curve
#' 
#' @description Calculates Saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' @param Tair     Air temperature (deg C)
#' @param formula  Formula to be used. Either "Alduchov_1996" or "Sonntag_1990"
#' 
#' @details Esat (kPa) is calculated based on the Magnus equation:
#' 
#'  \deqn{Esat = a * exp((b * Tair) / (c + Tair)) / 1000}
#'  
#'  where the coefficients a, b, c take different values depending on the formula used.
#'  The slope of the Esat curve (\eqn{\delta}) is calculated as the first derivative of the function:
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
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
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
#' @return the psychrometric constant \eqn{\gamma} (kPa K-1)
#'  
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
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
#' @return Latent heat of vaporization \eqn{\lambda} (J kg-1) 
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
#' @aliases LE.to.ET ET.to.LE
#' 
#' @param LE   Latent heat flux (W m-2)
#' @param ET   Evapotranspiration (kg m-2 s-1)
#' @param Tair Air temperature (deg C)
#' 
#' @export
LE.to.ET <- function(LE,Tair){

  Tair   <- Tair
  lambda <- LE.vaporization(Tair)
  ET     <- LE/lambda
   
  return(ET)
}

ET.to.LE <- function(ET,Tair){
  
  Tair   <- Tair
  lambda <- LE.vaporization(Tair)
  LE     <- ET*lambda
  
  return(LE)
}


#' Conversion between conductance units
#' 
#' @description Converts conductances from ms-1
#'              to mol m-2 s-1, or vice versa.
#' 
#' @aliases ms.to.mol mol.to.ms
#' 
#' @param G_ms       Conductance (ms-1)
#' @param G_mol      Conductance (mol m-2 s-1)
#' @param Tair       Air temperature (deg C)
#' @param pressure   Air pressure (kPa)
#' @param constants  Kelvin - conversion degree Celsius to Kelvin \cr
#'                   Rgas - universal gas constant (J mol-1 K-1)
#' 
#' @family conductance conversion
#' 
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology.
#'             2nd Edition., 2nd Edn. Cambridge University Press, Cambridge. 428 p --> replace with third edition
#'             
#' @export
ms.to.mol <- function(G_ms,Tair,pressure,constants=bigleaf.constants()){
  Tair   <- Tair + constants$Kelvin
  G_mol  <- G_ms * pressure*1000 / (constants$Rgas * Tair)
  
  return(G_mol)
}


#' @rdname ms.to.mol
#' @family conductance conversion
mol.to.ms <- function(G_mol,Tair,pressure,constants=bigleaf.constants()){
  Tair  <- Tair + constants$Kelvin
  G_ms  <- G_mol * (constants$Rgas * Tair) / (pressure*1000)
  
  return(G_ms)
}



#' Conversions between humidity measures
#' 
#' @description Conversion between vapor pressure, vapor pressure deficit (VPD), specific humidity, and relative humidity.
#' 
#' @param e         Vapor pressure (kPa)
#' @param p         Air pressure (kPa)
#' @param q         Specific humidity (kg kg-1)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param rH        Relative humidity (-)
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





#' Conversions between radiation measures
#'
#'
#'
mol.to.Watt <- function(PPFD,conversion){
  PAR <- PPFD * conversion
  return(PAR)
}

Watt.to.mol <- function(PAR,conversion){
  PPFD <- PAR
  return(PPFD)
}


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



###########################
#### Evapotranspiration ###-------------------------------------------------------------------------------
###########################

#' Potential evapotranspiration from the Priestley-Taylor equation
#' 
#' @details 
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2)
#' @param S         Sum of all storage fluxes (W m-2)
#' @param alpha     Priestley-Taylor coefficient (-)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @details Potential evapotranspiration is calculated according to Priestley & Taylor, 1972:
#' 
#'          \deqn{(\alpha * \delta * (Rn - G - S)) / (\delta + \gamma)}
#'
#' @return a data.frame with the following columns:
#'         \item{ET_pot}{potential evapotranspiration (kg m-2 s-1)}
#'         \item{LE_pot}{potential latent heat flux (W m-2)}
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
  ET_pot <- LE.to.ET(LE_pot,Tair)
  
  return(data.frame(ET_pot,LE_pot))
}



#' Imposed and Equilibrium Evapotranspiration
#' 
#' @description Evapotranspiration (ET) split up in imposed ET and equilibrium ET.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param VPD       Air vapor pressure deficit (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param Gs        surface conductance (m s-1)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' 
#' @details Total evapotranspiration can be written in the the form:
#' 
#'          \deqn{ET = \Omega ET_eq + (1 - \Omega)ET_imp}
#'          
#'          where ET_eq is the equilibrium evapotranspiration rate, the ET rate 
#'          that would occur under uncoupled conditions, where the heat budget
#'          is dominated by radiation.
#'          
#'          \deqn{ET_eq = \delta * (Rn - G - S) * \lambda / (\delta + \gamma)}
#'          
#'          and ET_imp is the imposed evapotranspiration rate, the ET rate
#'          that would occur under fully coupled conditions:
#'          
#'          \deqn{ET_imp = \rho * cp * VPD * Gs * \lambda / \gamma}
#' 
#' @note Surface conductance (Gs) is calculated with \code{\link{surface.conductance}}.
#'       Aerodynamic conductance (Ga) can be calculated using \code{\link{Ga}}.
#'       
#' @return a data.frame with the following columns:
#'         \item{ET_eq}{equilibrium ET (kg m-2 s-1)}
#'         \item{ET_imp}{imposed ET (kg m-2 s-1)}
#'         \item{LE_eq}{equilibrium LE (W m-2)}
#'         \item{LE_imp}{imposed LE (W m-2)}      
#' 
#' @references Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49.
#'             
#'             Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#' @return
ET.components <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",
                          Rn="Rn",Gs="Gs",constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  VPD      <- check.columns(data,VPD)
  Rn       <- check.columns(data,Rn)
  Gs       <- check.columns(data,Gs)
  
  rho    <- air.density(Tair,pressure)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  delta  <- Esat(Tair)[,"delta"]
  
  LE_eq  <- (delta * Rn) / (gamma + delta)
  LE_imp <- (rho * constants$cp * Gs * VPD) / gamma
  
  ET_imp <- LE.to.ET(LE_imp,Tair)
  ET_eq  <- LE.to.ET(LE_eq,Tair)
  
  return(data.frame(ET_eq,ET_imp,LE_eq,LE_imp))
}


#' Bulk intercellular CO2 concentration
#' 
#' @description Bulk canopy intercellular CO2 concentration (Ci) calculated based on Fick's law
#'              given surface conductance (Gs), gross primary productivity (GPP) and 
#'              atmospheric CO2 concentration (Ca).
#'              
#' @param Ca         Atmospheric CO2 concentration (umol mol-1)              
#' @param GPP        Gross primary productivity (umol CO2 m-2 s-1)              
#' @param Gs         Surface conductance to water (mol m-2 s-1)
#' @param surface.Ca Should the derived surface CO2 concentration be used instead of 
#'                   measured atmospheric CO2? If TRUE, Ca is derived as shown in \code{Details}.
#' @param Ga         Aerodynamic conductance to CO2 (mol m-2 s-1) 
#' @param NEE        Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param constants  DwDc Ratio of the molecular diffusivities for water vapor and CO2 (-)
#' 
#' @details Bulk intercellular CO2 concentration (Ci) is given as:
#' 
#'          \deqn{Ci = Ca - GPP/(Gs/1.6)}
#'          
#'          where Gs/1.6 represents the surface conductance to CO2.
#'          
#' @note The equation is based on Fick's law of diffusion and is equivalent to the
#'       often used equation at leaf level (ci = ca - An/gs).
#'       Note that GPP and Gs have a different interpretation than An and gs.
#'       Gs comprises non-pyhsiological contributions (i.e. physical evaporation)
#'       and is confounded by physical factors (e.g. energy balance non-closure).
#'       GPP does not account for dark respiration and is further subject to uncertainties
#'       in the NEE partitioning algorithm used.
#'       This function should be used with care and the resulting Ci might not be
#'       readily comparable to its leaf-level analogue.          
#' 
#' @return Ci - Bulk canopy intercellular CO2 concentration (umol mol-1)
#' 
#' 
#' @export
intercellular.CO2 <- function(Ca,GPP,Gs,surface.Ca=F,Ga=NULL,NEE=NULL,
                              constants=bigleaf.constants()){
  
  if (surface.Ca){
    Ca <- Ca.surface(Ca,NEE,Ga)
  }
  
  Ci <- Ca - GPP/(Gs/constants$DwDc)
  
  return(Ci)
  
}







########################
#### Energy balance ####-------------------------------------------------------------------------------
########################

#' Biochemical energy
#' 
#' @description Radiant energy absorbed in photosynthesis or heat release by respiration calculated
#'              from net ecosystem exchange of CO2 (NEE).  
#' 
#' @param alpha   Energy taken up/released by photosynthesis/respiration (J umol-1)
#' @param NEE     Net ecosystem exchange (umol CO2 m-2 s-1)
#' 
#' @details The following sign convention is employed: NEE is negative when carbon is taken up by the ecosystem.
#'          Positive values of the resulting biochemical energy mean that energy is taken up by the ecosystem, 
#'          negative ones that heat is released.
#'          The value of alpha is taken from Nobel 1974 (see Meyers & Hollinger 2004), but other values
#'          have been used (e.g. Blanken et al., 1997)
#' 
#' @return Sp - biochemical energy (W m-2)
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
#' @param data Data.frame or matrix containing all required variables
#' @param Rn   Net radiation (W m-2)
#' @param G    Ground heat flux (W m-2); optional
#' @param S    Sum of all storage fluxes (W m-2); optional
#' @param LE   Latent heat flux (W m-2)
#' @param H    Sensible heat flux (W m-2)
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

