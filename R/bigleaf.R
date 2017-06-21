#### Functions for R-package 'bigleaf' ------------------------------------------------------------


# library(devtools)
# library(roxygen2)


## Todo
# - how to define "families"
# - evtl new function for quality control
# - function for simple Gs calculation

# Notes:
# other estimates for d and z0m: Shaw and Pereira 1982 (see also Yang_2003)
# v = kinematic viscosity of the air (m^2/s); from=Massman_1999b ; compare with http://www.engineeringtoolbox.com/dry-air-properties-d_973.html


# ## artificial test data
# wind     <- rnorm(20,3,0.8)
# ustar    <- rnorm(20,0.5,0.2)
# Tair     <- rnorm(20,25,2)
# pressure <- rep(101.325,20)
# Rn       <- rnorm(20,400,40)
# H        <- 0.4*Rn
# LE       <- 0.5*Rn
# G        <- 0.1*Rn
# Ga       <- rnorm(20,0.08,0.01)
# Gs       <- rnorm(20,0.3,0.02)*(1/41)
# VPD      <- rnorm(20,2,0.3)
# precip   <- c(rep(0,18),0.1,0.001)
# PPFD     <- seq(180,300,length.out=20)
# GPP      <- rnorm(20,20,3)
# NEE      <- -GPP
# day      <- c(rep(1,10),rep(2,10))
# #
# #
# #
# testdf <- data.frame(wind,ustar,Tair,pressure,Rn,
#                      H,LE,G,Ga,Gs,VPD,precip,PPFD,GPP,NEE,day)


## real data
# path.pkg <- "C:/Profiles/jknauer/Documents/RPackage_BigLeaf/"
# load(paste0(path.pkg,"DK-Sor_2010_processed.rda"))
# testdk <- DK_Sor_2010_processed
# colnames(testdk)[10] <- "VPD"
# testdk[,"VPD"] <- testdk[,"VPD"]/10
# save(testdk,file="C:/Profiles/jknauer/Desktop/testdk.RData")


########################
### helper functions ###
########################

#' checks columns in a data frame 
#' 
#' @description check if columns in a data.frame or matrix (or vector variables) 
#'              exist and if they are of the right type (numeric).
#' 
#' @param data         a data.frame or a matrix
#' @param column_name  the column name
#' 
#' 
#' 
check.columns <- function(data,column_name){

  if (is.character(column_name)){
    if (!missing(data)){
      if (length(column_name) == 1){
        if (column_name %in% colnames(data)){
          var <- data[,column_name]
          if (is.numeric(var)){
            return(unname(var))
          } else {
            stop("variable '",column_name,"' must be numeric",call.=FALSE)
          }
        } else {
          stop ("variable '",column_name,"' does not exist in the input matrix/data.frame",call.=FALSE)
        }
      } else {
        stop("name of variable '",deparse(substitute(column_name)),"' must have length 1",call.=FALSE)
      } 
    } else {
      if ("data" %in% names(formals(sys.function(which=-1)))){
        stop("variable '",deparse(substitute(column_name)),"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",column_name,"' as a numeric vector, or an input matrix/data.frame with a column named '",column_name,"'",call.=FALSE)
      } else {
        stop("variable '",deparse(substitute(column_name)),"' must be numeric",call.=FALSE)
      }
    }
  } else {
    if (!missing(data)){
      if (is.numeric(column_name) & length(column_name) == nrow(data)){
        return(unname(column_name))
      } else if (is.numeric(column_name) & length(column_name) != nrow(data)) {
        stop("variable '",deparse(substitute(column_name)),"' must have the same number of observations as the input matrix/data.frame",call.=FALSE)
      } else if (!is.numeric(column_name)){
        stop("variable '",column_name,"' must be numeric",call.=FALSE)
      }
    } else {
      if (is.numeric(column_name)){
        return(unname(column_name))
      } else {
        stop("variable '",column_name,"' must be numeric",call.=FALSE)
      } 
    } 
  }
}


#' tests if variables have the same length
#' 
#' 
check.length <- function(...){

  args <- list(...)
  length.args <- sapply(args,length)
  
  if (length(unique(length.args)) > 1){
    stop("all input variables must have the same length",call.=FALSE)
  }
}


  
#########################
### global constants ####
#########################

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


########################
### Filter functions ###----------------------------------------------------------------------------
########################

#' Filter data
#'
#' @description Filters timeseries in which Gs constitutes a physiologically meaningful
#'              quantity.
#' 
#' @param data          data.frame or matrix containing all required input variables in 
#'                      half-hourly or hourly resolution. Including year, month, day information
#' @param precip        precipitation (mm)
#' @param PPFD          photosynthetically photon flux density (umol m-2 s-1)
#' @param Tair          air temperature (deg C)
#' @param ustar         friction velocity (m s-1)
#' @param VPD           vapor pressure deficit (kPa)
#' @param GPP           gross primary productivity (umol m-2 s-1)
#' @param doy           day of year
#' @param year          year
#' @param quality.control should quality control be applied?
#' @param quality.ext   the extension to the variables' names that marks them as quality control variables
#' @param vars.qc       character vector indicating the variables for which quality filter should 
#'                      be applied. Leave empty if no quality filter should be applied.  
#' @param good.quality  which values indicate good quality in the quality control (qc) variables?
#' @param missing.qc.as.bad If quality control variable is NA, should the corresponding data point be
#'                          treated as bad quality? Defaults to TRUE
#' @param tprecip       precipitation threshold used to demark a precipitation event (mm)
#' @param precip.hours  number of hours removed following a precipitation event (h)
#' @param NA.as.precip  if TRUE missing precipitation data are treated as precipitation events
#' @param trad          radiation threshold (umol m-2 s-1)
#' @param ttemp         Temperature threshold (deg C)
#' @param tustar        Friction velocity threshold (m s-1)
#' @param tGPP          GPP threshold (fraction of 95% quantile of the GPP time series).
#'                      Must be between 0 and 1. 
#' @param ws            Window size used for GPP time series smoothing
#' @param min.int       Minimum time interval in days for a given state of growing season
#' @param trH           Relative humidity threshold (-). Note that relative humidity
#'                      is calculated from VPD internally

#' 
#' 
#' @details This routine consists of two parts:
#'          1) Quality control: All variables included in \code{vars.qc} are filtered for 
#'             good quality data. For these variables a corresponding quality variable with 
#'             the same name as the variable plus the extension as specified in \code{quality.ext}
#'             must be provided. For timesteps where the value of the quality indicator is not included
#'             in the argument \code{good.quality}, i.e. the quality is not considered as 'good', 
#'             its value is set to NA.
#'             
#'          2) Meteorological filtering. Under certain conditions (e.g. low ustar), the assumptions
#'             of the EC method are not fulfilled. Further, some data analysis require certain meteorological
#'             conditions, such as precipitation free periods, or active vegetation (growing season, daytime).
#'             The filter applied in this second step serves to exclude time periods that do not fulfill the criteria
#'             specified in the arguments. More specifically, timeperiods where one of the variables is lower than the 
#'             specified threshold are set to NA for all variables. E.g. when radiation is below 'trad', all variables 
#'             are set to NA.
#'             Set a threshold to 0 (or the minimum of the variable) if the dataset should not be filtered for a specific variable.
#'          
#' @return filtered_data The same data frame as provided as first argument to the function
#'                       with an additional column "valid", which indicates whether the timestep is valid 
#'                       according to the filter criteria (1) or invalid (0).
#' 
#' @export                     
filter.data <- function(data,precip="precip",PPFD="PPFD",Tair="Tair",ustar="ustar",
                        VPD="VPD",GPP="GPP_nt",doy="doy",year="year",quality.control=TRUE,
                        quality.ext="_qc",vars.qc=c("precip","Tair","VPD","GPP_nt","H","LE"),
                        good.quality=c(0,1),missing.qc.as.bad=TRUE,
                        tprecip=0.01,precip.hours=24,NA.as.precip=F,trad=200,ttemp=5,tustar=0.2,tGPP=0.5,
                        ws=15,min.int=5,trH=0.95){
  
  
  ### I) Quality control filter
  if (quality.control){
    
    if (any(!vars.qc %in% colnames(data))){
      
      missing_vars <- vars.qc[which(!vars.qc %in% colnames(data))]
      stop(paste("Variable ",missing_vars," is included in 'vars.qc', but does not exist in the input data!"))
      
    }
    
    vars.qc_qc <- paste0(vars.qc,quality.ext)
    if (any(!vars.qc_qc %in% colnames(data))){
      
      missing_vars_qc <- vars.qc_qc[which(!vars.qc_qc %in% colnames(data))]
      missing_vars2   <- substr(missing_vars_qc,1,nchar(missing_vars_qc) - nchar(quality.ext))
      stop(paste("Quality control for variable ",missing_vars2,"(",missing_vars_qc,") does not exist in the input data!")) 
    }

    ## data quality
    cat("Quality control:",fill=TRUE)
    for (var in vars.qc){
      assign(var,check.columns(data,var))
      assign(paste0(var,quality.ext),check.columns(data,paste0(var,quality.ext))) # create variable "*_qc", make quality check before
      
      if (missing.qc.as.bad){
        data[get(paste0(var,quality.ext)) > max(good.quality) | is.na(get(paste0(var,quality.ext))),var] <- NA   # exclude bad quality data or those where qc flag is not available
        qc_invalid      <- sum(get(paste0(var,quality.ext)) > max(good.quality) | is.na(get(paste0(var,quality.ext)))) # count & report
      } else { # same, but consider missing quality flag variables as good
        data[get(paste0(var,quality.ext)) > max(good.quality) & !is.na(get(paste0(var,quality.ext))),var] <- NA
        qc_invalid      <- sum(get(paste0(var,quality.ext)) > max(good.quality) & !is.na(get(paste0(var,quality.ext))))
      }
      
      qc_invalid_perc <- round((qc_invalid/nrow(data))*100,2)
        
      cat(var,": ",qc_invalid," data points (",qc_invalid_perc,"%) set to NA",fill=T,sep="")
    }
    cat("-------------------------------------",fill=T)
  }
  
  ## test data availability
  precip <- check.columns(data,precip)
  PPFD   <- check.columns(data,PPFD)
  Tair   <- check.columns(data,Tair)
  ustar  <- check.columns(data,ustar)
  VPD    <- check.columns(data,VPD)
  GPP    <- check.columns(data,GPP)
  doy    <- check.columns(data,doy)
  year   <- check.columns(data,year)
  
  date <- strptime(paste0(year,"-",doy),format="%Y-%j")
  
  #### II) filtering based on meteorology
  rH <- VPD.to.rH(VPD,Tair)
  
  ## hourly or halfhourly data?
  tstep_day <- as.integer(nrow(data) / length(unique(year)))
  thour     <- ifelse(tstep_day >= 17520 & any(!is.na(Tair[seq(1,length(Tair),2)])) &
                      any(!is.na(Tair[seq(2,length(Tair),2)])),2,1)
  
  valid <- rep(1L,nrow(data))
  
  
  ## start filtering
  # 1) GPP
  GPP_daily        <- aggregate(GPP,by=list(strftime(date)),mean,na.rm=T)[,2]
  growing_season   <- filter.growseas(GPP_daily,tGPP=tGPP,ws=ws,min.int=min.int)
  growseas_invalid <- which(sapply(growing_season,rep,48) == 0)
  
  # 2) precipitation
  if (NA.as.precip){
    precip_events <- which(precip > tprecip | is.na(precip))
  } else {
    precip_events <- which(precip > tprecip)
  }
  precip_invalid <- unique(as.numeric(unlist(sapply(precip_events, function(x) x:(min(x+precip.hours*thour,nrow(data),na.rm=T))))))
  
  # 3) meteorological variables (PPFD, Tair, ustar, rH) 
  PPFD_invalid  <- which(PPFD <= trad | is.na(trad))
  Tair_invalid  <- which(Tair <= ttemp | is.na(Tair))
  ustar_invalid <- which(ustar <= tustar | is.na(ustar))
  rH_invalid    <- which(rH >= trH | is.na(rH))
  
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
  
  return(data_filtered)
}




#' GPP-based growing season filter
#' 
#' @description Filters annual time series for growing season based on smoothed daily GPP data.
#' 
#' @param GPPd    daily GPP (any unit) 
#' @param tGPP    GPP threshold (fraction of 95% quantile of the GPP time series).
#'                Takes values between 0 and 1. 
#' @param ws      window size used for GPP time series smoothing
#' @param min.int minimum time interval in days for a given state of growing season
#' 
#' @details The basic idea behind the growing season filter is that vegetation is 
#'          considered to be active when its carbon uptake (GPP) is above a specified 
#'          threshold, which is defined relative to the peak GPP (95th-percentile) 
#'          observed in the year. 
#'          The GPP-threshold is calculated as:
#'          
#'          \deqn{GPP_threshold = quantile(GPPd,0.95)*tGPP}
#'          
#'          GPPd time series are smoothed with a moving average to avoid fluctuations 
#'          in the delineation of the growing season. The window size defaults to 15 
#'          days, but depending on the ecosystem, other values can be appropriate. 
#'          
#'          The argument \code{min.int} serves to avoid short fluctuations in the 
#'          status growing season vs. no growing season by defining a minimum length
#'          of the status. If a time interval shorter than \code{min.int} is labelled
#'          as growing season or non-growing season, it is changed to the status of 
#'          the neighbouring values.
#'          
#' @return a vector of type integer of the same length as the input GPPd in which 0 indicate
#'         no growing season (dormant season) and 1 indicate growing season.
#'                         
#' @export  
filter.growseas <- function(GPPd,tGPP,ws=15,min.int=5){
  
  if(sum(is.na(GPPd)) < 0.5*length(GPPd)){
  
    growseas      <- rep(1,length(GPPd))
    GPP_threshold <- quantile(GPPd,probs=0.95,na.rm=TRUE)*tGPP
    
    ## smooth GPP
    GPPd_smoothed <- filter(GPPd,method="convolution",filter=rep(1/ws,ws))
    
    ## set values at the beginning and end of the time series to the mean of the original values
    wsd <- floor(ws/2)
    GPPd_smoothed[1:wsd] <- mean(GPPd[1:(2*wsd)],na.rm=T)
    GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] <- mean(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)],na.rm=T)
    
    # check for occurence of missing values and set them to mean of the values surrounding them
    missing <- which(is.na(GPPd_smoothed))
    if (length(missing) > 0){
      if (length(missing) > 10){warning("Attention, there is a gap in 'GPPd' of length n = ",length(missing))}
      replace_val <- mean(GPPd_smoothed[max(1,missing[1] - 4):min((missing[length(missing)] + 4),length(GPPd_smoothed))],na.rm=T)
      GPPd_smoothed[missing] <- replace_val
    }
    
    # filter daily GPP
    growseas[GPPd_smoothed < GPP_threshold] <- 0
    
    ## change short intervals to the surrounding values to avoid 'wrong' fluctuations
    intervals <- rle(growseas)
    short_int <- which(intervals$lengths <= min.int)
    
    if (length(short_int) > 0){
      start <- numeric()
      end   <- numeric()
      
      for (i in 1:length(short_int)){
        
        start[i] <- sum(intervals$lengths[1:short_int[i]-1]) + 1
        end[i]   <- start[i]+intervals$lengths[short_int[i]] - 1

        val <- unique(growseas[start[i]:end[i]])
        
        if (val == 0 & growseas[start[i]-1] == 1){
          growseas[start[i]:end[i]] <- 1   
        } else if (val == 1 & growseas[start[i]-1] == 0){
          growseas[start[i]:end[i]] <- 0
        }
      }
    }

    growseas <- as.integer(growseas)
  
  } else {
    
    warning("number of available GPPd data is less than half the total number of days per year. Filter is not applied!")
    growseas <- as.integer(rep(1,length(GPPd)))
  
  }
  
  return(growseas)
}





###############################
### WUE metrics calculation ###
###############################

#' Water-Use Efficiency metrics
#' 
#' @description Calculation of various water use efficiency (WUE) metrics.
#' 
#' @param data      data.frame or matrix containing all required variables
#' @param GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' @param NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Tair      Air temperature (deg C)
#' @param constants Cmol - molar mass of carbon (kg mol-1)
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
#'          All metrics are calculated based on the median of all values. E.g.
#'          WUE = median(GPP/ET,na.rm=T)
#' 
#' @return a named vector with the following elements:
#'         \item{WUE}{Water-use efficiency (gC (kg H20)-1)}
#'         \item{WUE_NEE}{Water-use efficiency based on NEE (gC (kg H20)-1)}
#'         \item{IWUE}{Inherent water-use efficiency (gC kPa (kg H20)-1)}
#'         \item{uWUE}{Underlying water-use efficiency (gC kPa (kg H20)-1)}
#' 
#' @note Units for VPD can also be hPa. Units change accordingly.
#'       WUE_NEE is calculated based on the absolute value of NEE (the sign convention does not matter).
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
  check.length(GPP,NEE,LE,VPD,Tair)
  
  ET  <- LE.to.ET(LE,Tair)                 # kg H2O s-1
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




################################################
### Boundary layer conductance formulations #### 
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
#'  \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'  The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
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
#' @examples 
#' # calculate Canopy boundary layer resistance when ustar=0.4ms-1
#' Gb.Thom(0.4)
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
#' @details The Roughness Reynolds Number is calculated as in Massman 1999a:
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
#' @param Tair      Air temperature (degC)
#' @param presssure Atmospheric pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param wind      Wind speed (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zh        Canopy height (m)
#' @param zr        Reference height (m)
#' @param d         Zero-plane displacement height (-), can be calculated using \code{roughness.parameters}
#' @param Dl        Leaf characteristic dimension (m)
#' @param fc        Fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI       One-sided leaf area index (-)
#' @param N         Number of leaf sides participating in heat exchange (defaults to 2)
#' @param Cd        Foliage drag coefficient (-)
#' @param hs        Roughness height of the soil (m)
#' @param stab_formulation Stability correction function used (If stab_correction is TRUE).
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
#'          \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'          The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
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
#' @examples 
#' # Canopy boundary layer resistance (and kB-1 parameter) for a set of meteorological conditions,
#' # a leaf characteristic dimension of 1cm, and an LAI of 5
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
#' Gb.Su(data=df,zh=25,zr=40,Dl=0.01,LAI=5)
#' 
#' # the same meteorological conditions, but larger leaves
#' Gb.Su(data=df,zh=25,zr=40,d=17.5,Dl=0.1,LAI=5)
#' 
#' # same conditions, large leaves, and sparse canopy cover (LAI = 1.5)
#' Gb.Su(data=df,zh=25,zr=40,d=17.5,Dl=0.01,LAI=1.5)
#' @export
Gb.Su <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",wind="wind",
                  H="H",zh,zr,d,Dl,fc=NULL,LAI=NULL,N=2,Cd=0.2,hs=0.01,
                  stab_formulation=c("Dyer_1970","Businger_1971"),
                  constants=bigleaf.constants()){
        
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  wind     <- check.columns(data,wind)
  H        <- check.columns(data,H)
  check.length(Tair,pressure,ustar,wind,H)
  
  if (is.null(fc)) {
    if (is.null(LAI)){
      stop("one of 'fc' or 'LAI' must be provided",call.=FALSE)
    } else {
      fc  <- (1-exp(-LAI/2)) 
    }
  } 
  
  wind_zh <- wind.profile(data=data,heights=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,zh=zh,d=d,stab_correction=TRUE,
                          stab_formulation=stab_formulation)[,1]

  v   <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"v"]
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)[,"Re"]
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind_zh / v
  Ct  <- 1*0.71^-0.6667*Reh^-0.5*N   # 0.71 = Prandtl number

  kB     <- (constants$k*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb     <- kB/(constants$k*ustar)
  Rb_CO2 <- constants$Rbwc * Rb
  Gb     <- 1/Rb


  return(data.frame(Rb,Rb_CO2,Gb,kB))
}






#' Boundary layer conductance according to Choudhury & Monteith 1988
#' 
#' @param wind             Wind speed at sensor height (m s-1)
#' @param ustar            Friction velocity (m s-1)
#' @param leafwidth        Leaf width (m)
#' @param LAI              One-sided leaf area index
#' @param zh               Canopy height (m)
#' @param zr               Instrument (reference) height (m)
#' @param d                Zero-plane displacement height (-), can be calculated using \code{roughness.parameters}
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
#'          \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'          The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
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
#' @examples 
#' ## bulk canopy boundary layer resistance for a closed canopy (LAI=5) with large leaves (leafwdith=0.1)            
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250))    
#' Gb.Choudhury(data=df,leafwidth=0.1,LAI=5,zh=25,d=17.5,zr=40)
#' 
#' ## same conditions, but smaller leaves (leafwidth=0.01)
#' Gb.Choudhury(data=df,leafwidth=0.01,LAI=5,zh=25,d=17.5,zr=40) 
#' 
#' @export                                                                                                                                                                                                                                                                                    
Gb.Choudhury <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                         leafwidth,LAI,zh,zr,d,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  wind     <- check.columns(data,wind)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)

  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)
  wind_zh <- wind.profile(data=data,heights=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,zh=zh,d=d,stab_correction=TRUE,
                          stab_formulation=stab_formulation)[,1]
  
  ## avoid zero windspeed
  wind_zh <- pmax(0.01,wind_zh)
  
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
#' @param data      Data.frame or matrix containing all required variables
#' @param heights   vector with heights for which wind speed is to be 
#'                  calculated.
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)                                                                                  
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zr        Instrument (reference) height (m)
#' @param zh        Canopy height (m)
#' @param d         Zero-plane displacement height (-)
#' @param frac_d    Fraction of displacement height on canopy height (-);
#'                  only used if \code{d} is not available
#' @param z0m       Roughness length (m), optional; only used if stab_correction=FALSE (default=0.1) 
#' @param frac_z0m  Fraction of roughness length on canopy height (-), optional; only used 
#'                  if stab_correction=FALSE (default=0.1), only used if \code{z0m} is not available
#' @param stab_correction Should stability correction be applied? Defaults to TRUE
#' @param stab_formulation Stability correction function used (If stab_correction is TRUE).
#'                         Either "Dyer_1970" or "Businger_1971".
#' @param constants k - von-Karman constant (-) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  g - gravitational acceleration (m s-2) \cr
#'                                   
#' @note Note that this equation is only valid for z >= d + z0m, and it is not 
#'       meaningful to calculate values closely above d + z0m. All values in \code{heights}
#'       smaller than d + z0m will return 0.                                 
#'                                  
#' @return a data.frame with rows representing time and columns representing heights as specified in \code{heights}.     
#'                                            
#' @references
#' 
#' @examples 
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
#' wind.profile(df,heights=c(17:40),zr=40,zh=25,d=16)
#' 
#' @export                                                                                                                          
wind.profile <- function(data,heights,Tair="Tair",pressure="pressure",ustar="ustar",
                         H="H",zr,zh,d=NULL,frac_d=0.7,z0m=NULL,frac_z0m=0.1,stab_correction=TRUE,
                         stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  check.length(Tair,pressure,ustar,H)
  
  wind_heights <- data.frame(matrix(NA,ncol=length(heights),nrow=length(Tair)))
  colnames(wind_heights) <- paste0(heights,"m")
  
  if (stab_correction){
    
    if (is.null(d)){
      d <- frac_d * zh
    }
        
    z0m <- roughness.parameters(method="wind_profile",zh=zh,zr=zr,d=d,data=data,
                                stab_roughness=TRUE,constants=constants)[,"z0m"]
    
    if ( any(heights < (d + z0m) & !is.na(d + z0m)) ){
      warning("function is only valid for heights above d + z0m! Wind speed for heights below d + z0m will return 0!") 
    } 
  }  
    
  for (z in heights){
    i <- which(heights == z)

    if (stab_correction){
      
      zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=z,d=d,constants=constants)
      
      psi_m <- stability.correction(zeta)[,"psi_m"]
      
      wind_heights[,i] <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m) - psi_m))
    
    } else {
      
      if (is.null(d)){
        d <- frac_d * zh
      }
      
      if (is.null(z0m)){
        z0m <- frac_z0m * zh
      }

      
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
#' @param method    Method to use, either "canopy_height", "canopy_height&LAI", or "wind_profile" \cr
#'                  NOTE: if method is "canopy_height", only the following three arguments
#'                  are used. If method is "canopy_height&LAI", only zh, LAI, cd, and hs are required     
#' @param zh        Vegetation height (m)          
#' @param frac_d    Fraction of displacement height on canopy height (-)
#' @param frac_z0m  Fraction of roughness length on canopy height (-)
#' @param LAI       Leaf area index (-) 
#'                  The following arguments are only needed if \code{method = "wind_profile"}!
#' @param zr        Instrument (reference) height (m)
#' @param cd        Mean drag coefficient for individual leaves. Defaults to 0.2.
#' @param hs        roughness length of the soil surface (m)
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param wind      Wind speed at height zr (m s-1)
#' @param ustar     Friction velocity (m s-1)
#' @param H              Sensible heat flux (W m-2)
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
#'          Alternatively, d and z0m can be estimated from both canopy height and LAI
#'          (If \code{method} is \code{canopy_height&LAI}).
#'          Based on data from Shaw & Pereira 1982, Choudhury & Monteith 1988 proposed 
#'          the following semi-empirical relations:
#'          
#'          X = cd * L
#'          d = 1.1 * zh * ln(1 + X^(1/4)) 
#'          z0m = hs + 0.3 * zh * X^(1/2)  for 0 <= X <= 0.2
#'          z0m = hs * zh * (1 - d/zh)     for 0.2 < X 
#'          
#'          If \code{method} is \code{wind_profile}, z0m is estimated by solving
#'          the windspeed profile for z0m:
#'          
#'          \deqn{z0m = median((zr - d) * exp(-k*wind / ustar - psi_m)}
#'                  
#'          By default, d in this equation is fixed to 0.7zh, but can be set to any
#'          other value. psi_m is 0 if \code{stab_roughness=FALSE} (the default).       
#' 
#' @return a data.frame with the following columns:
#'         \item{d}{Zero-plane displacement height (m)}
#'         \item{z0m}{Roughness length for momentum (m)}
#'         \item{z0m_se}{Only if \code{method} is \code{wind_profile}: Standard Error of the median for z0m (m)}
#'
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             Shaw, R. H., Pereira, A., 1982: Aerodynamic roughness of a plant canopy: 
#'             a numerical experiment. Agricultural Meteorology, 26, 51-65.
#'    
#' @examples 
#' # estimate d and z0m from canopy height for a dense (LAI=5) and open (LAI=2) canopy
#' roughness.parameters(method="canopy_height&LAI",zh=25,LAI=5)
#' roughness.parameters(method="canopy_height&LAI",zh=25,LAI=2)   
#'    
#' # fix d to 0.7*zh and estimate z0m from the wind profile
#' df <- data.frame(Tair=c(25,25,25),pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=200)
#' roughness.parameters(method="wind_profile",zh=25,zr=40,frac_d=0.7,data=df)
#' 
#' # assume d = 0.8*zh
#' roughness.parameters(method="wind_profile",zh=25,zr=40,frac_d=0.8,data=df)  
#' @export                                  
roughness.parameters <- function(method=c("canopy_height","canopy_height&LAI","wind_profile"),zh,
                                 frac_d=0.7,frac_z0m=0.1,LAI,zr,cd=0.2,hs=0.01,data,Tair="Tair",pressure="pressure",
                                 wind="wind",ustar="ustar",H="H",d=NULL,z0m=NULL,
                                 stab_roughness=FALSE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf.constants()){
  
  method           <- match.arg(method)
  stab_formulation <- match.arg(stab_formulation)
  
  if (method == "canopy_height"){
  
    d      <- frac_d*zh
    z0m    <- frac_z0m*zh
    z0m_se <- NA
  
  } else if (method == "canopy_height&LAI"){
    
    X <- cd * LAI
    d <- 1.1 * zh * log(1 + X^(1/4))
    
    if (X >= 0 & X <= 0.2){
      z0m <- hs + 0.3 * X^(1/2)
    } else {
      z0m <- 0.3 * zh * (1 - d/zh)
    }
    z0m_se <- NA
    
  } else if (method == "wind_profile"){
    
    Tair     <- check.columns(data,Tair)
    pressure <- check.columns(data,pressure)
    wind     <- check.columns(data,wind)
    ustar    <- check.columns(data,ustar)
    H        <- check.columns(data,H)
    check.length(Tair,pressure,wind,ustar,H)
    
    if (is.null(d)){
      
      d <- frac_d * zh
    
    }
    
    if (stab_roughness){
      
      zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=zr,d=d,constants=constants)
      psi_m <- stability.correction(zeta)[,"psi_m"]
      
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar - psi_m)
      
    } else {
      
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar)
      
    }
    
    z0m_all[z0m_all > zh] <- NA
    
    z0m    <- median(z0m_all,na.rm=T)
    z0m_se <- 1.253 * (sd(z0m_all,na.rm=T) / sqrt(length(z0m_all[complete.cases(z0m_all)])))

  }

  return(data.frame(d,z0m,z0m_se))
}


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
#' @examples 
#' # air density at 25degC and standard pressure (101.325kPa)
#' air.density(25,101.325)
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
    Temp <- virtual.temp(Tair,q,constants)  + constants$Kelvin
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
#' @param data      Data.frame or matrix containing all required variables
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
stability.parameter <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                H="H",zr,d,constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  ustar    <- check.columns(data,ustar)
  H        <- check.columns(data,H)
  check.length(Tair,pressure,ustar,H)
  
  MOL  <- MoninObukhov.length(Tair,pressure,ustar,H,constants)
  zeta <- (zr - d) / MOL
  
  return(zeta)
  
}


#' Integrated stability correction functions for heat and momentum
#' 
#' @description dimensionless stability functions needed to correct deviations
#'              from the exponential wind profile under non-neutral conditions.
#'              
#' @param zeta         Stability parameter zeta (-)
#' @param formulation  Formulation for the stability function. Either "Dyer_1970", or "Businger_1971"
#'
#' @details The functions give the integrated form of the universal functions. They
#'          depend on the value of the stability parameter \eqn{\zeta},
#'          which can be calculated from the function \code{\link{stability.parameter}}.
#'          The integration of the universal functions is:
#'          
#'          \deqn{\psi = -x * zeta} 
#'          
#'          for stable atmospheric conditions (\eqn{\zeta} >= 0), and
#'          
#'          \deqn{\psi = 2 * log( (1 + y) / 2) }
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
#' @examples 
#' zeta <- seq(-2,0.5,0.01)
#' stability.correction(zeta)
#' stability.correction(zeta,formulation="Businger_1971")                          
#'             
#' @export   
stability.correction <- function(zeta,formulation=c("Dyer_1970","Businger_1971")){
  
  formulation  <- match.arg(formulation)
  
  zeta <- check.columns(,zeta)
  
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





#' Aerodynamic conductance
#' 
#' @description Bulk aerodynamic conductance, including options for the boundary layer conductance
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
#' @param wind_profile      Should Ga for resistance be calculated based on the logarithmic wind profile equation? 
#'                          Defaults to FALSE.
#' @param stab_correction   Should stability correction be applied? Defaults to TRUE. Ignored if \code{wind_profile} is FALSE                         
#' @param stab_formulation  Stability correction function. Either "Dyer_1970" (default) or
#'                          "Businger_1971". Ignored if one of \code{wind_profile} and \code{stab_correction} is FALSE
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
#' Aerodynamic conductance for heat (Ga_h) is calculated as:
#' 
#'  Ga_h = 1 / (Ra_m + Rb)
#'  
#'  where Ra_m is the aerodynamic resistance for momentum and Rb the canopy boundary
#'  layer resistance ('excess resistance').
#'  
#'  The aerodynamic resistance for momentum Ra_m is given by:
#'  
#'  \deqn{Ra_m = u/ustar^2}
#'  
#'  Note that this formulation accounts for changes in atmospheric stability, and does not require an
#'  additional stability correction function. 
#'  
#'  An alternative method to calculate Ra_m is provided
#'  (calculated if \code{wind_profile} set to TRUE):
#'  
#'  \deqn{Ra_m = (ln((zr - d)/z0m) - psi_h) / (k ustar)}
#'  
#'  If the roughness parameters z0m and d are unknown, they can be estimated using
#'  \link{\code{roughness.parameters}}. The argument \code{stab_formulation} determines the stability correction function used 
#'  to account for the effect of atmospheric stability on Ra_m (Ra_m is lower for unstable
#'  and higher for stable stratification). Stratification is based on a stability parameter zeta (z-d/L),
#'  where z = reference height, d the zero-plane displacement height, and L the Monin-Obukhov length, 
#'  calculated with \code{\link{MoninObukhov.length}}
#'  The Stability correction function is chosen by the argument \code{stab_formulation}. Options are 
#'  Dyer_1970" and "Businger_1971".
#'  
#'  The model used to determine the canopy boundary layer resistance (Rb) is specified by 
#'  the argument \code{Rb_model}. The following options are implemented:
#'  "Thom_1972" is an empirical formulation based on the friction velocity (ustar) (Thom 1972):
#'  
#'    \deqn{Rb = 0.67ustar^0.67}
#'    
#'  The model by Choudhury 1988 ("Choudhury_1988"), calculates Rb
#'  based on leaf width, LAI and ustar (Note that in the original formulation leaf width
#'  instead of the characteristic leaf dimension (Dl) is taken):
#'   
#'     \deqn{Gb = LAI((2a/\alpha)*sqrt(u(h)/w)*(1-exp(-\alpha/2)))}
#'     
#'  The option "Su_2001" calculates Rb based on the physically-based Rb model by Su et al. 2001,
#'  a simplification of the model developed by Massmann 1999:
#'  
#'     \deqn{kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2}
#'  
#'  The models calculate the parameter kB-1, which is related to Rb:
#'  
#'     \deqn{kB-1 = Rb * (k * ustar)}
#'  
#'  Both kB-1 and Rb, as well as Gb are given in the output.
#'  
#'  Rb (and Gb) for water vapor and heat are assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'  \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'  The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
#'  It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'  across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'  mixing (Nobel 2005).
#'  
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
#' @note The roughness length for water and heat (z0h) is not returned by this function, but 
#'       it can be calculated from the following relationship (e.g. Verma 1989):
#'       
#'       kB-1 = ln(z0m/z0h) 
#'       
#'       it follows:
#'       
#'       z0h = z0m / exp(kB-1)    
#'       
#'       kB-1 is an output of this function
#'         
#' @references Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
#'             In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
#'             
#'             Verhoef, A., De Bruin, H., Van Den Hurk, B., 1997: Some practical notes on the parameter kB-1
#'             for sparse vegetation. Journal of Applied Meteorology, 36, 560-572.
#' 
#' @examples 
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250))   
#'  
#' # simple calculation of Ga  
#' aerodynamic.conductance(df,Rb_model="Thom_1972") 
#' 
#' # calculation of Ga using a model derived from the logarithmic wind profile
#' aerodynamic.conductance(df,Rb_model="Thom_1972",zr=40,zh=25,d=17.5,z0m=2,wind_profile=TRUE) 
#' 
#' # simple calculation of Ga, but a physically based canopy boundary layer model
#' aerodynamic.conductance(df,Rb_model="Su_2001",zr=40,zh=25,d=17.5,Dl=0.05,N=2,fc=0.8)
#' 
#' @export
aerodynamic.conductance <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                                    zr,zh,d,z0m,Dl,N,fc=NULL,LAI,Cd=0.2,hs=0.01,wind_profile=FALSE,
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
  check.length(Tair,pressure,wind,ustar,H)

  ## calculate boundary layer conductance
  if (Rb_model != "constant_kB"){
    
    if (Rb_model == "Thom_1972"){
      
      Gb_mod <- Gb.Thom(ustar=ustar,constants=constants)
      
    } else if (Rb_model == "Choudhury_1988"){
      
      Gb_mod <- Gb.Choudhury(data,leafwidth=Dl,LAI=LAI,zh=zh,zr=zr,d=d,
                             stab_formulation=stab_formulation,constants=constants)
      
    } else if (Rb_model == "Su_2001"){
      
      Gb_mod <- Gb.Su(data=data,zh=zh,zr=zr,d=d,Dl=Dl,N=N,fc=fc,LAI=LAI,
                      Cd=Cd,hs=hs,stab_formulation=stab_formulation,constants=constants)  
      
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
  
  if (wind_profile){
    
    if (stab_correction){
    
      zeta  <-  stability.parameter(data=data,zr=zr,d=d,constants=constants)
    
      if (stab_formulation %in% c("Dyer_1970","Businger_1971")){
      
        psi_h <- stability.correction(zeta)[,"psi_h"]
      
      } else {
        stop("'stab_formulation' has to be one of 'Dyer_1970' or 'Businger_1971'.
               Choose 'stab_correction = FALSE' if no stability correction should be applied.")
      }
  
      Ra_m  <- pmax((log((zr - d)/z0m) - psi_h),0) / (constants$k*ustar)
  
    } else {
    
      Ra_m  <- pmax((log((zr - d)/z0m)),0) / (constants$k*ustar)
    
      zeta = psi_h <- rep(NA,length=length(Ra_m))
    
    }

  } else {
    
    Ra_m <- wind / ustar^2
    
    zeta = psi_h <- rep(NA,length=length(Ra_m))
    
  }
  
  Ga_m   <- 1/Ra_m
  Ra_h   <- Ra_m + Rb
  Ra_CO2 <- Ra_m + Rb_CO2
  Ga_h   <- 1/Ra_h
  Ga_CO2 <- 1/Ra_CO2

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Ga_CO2,Ra_CO2,Rb,Rb_CO2,Gb,kB,zeta,psi_h))
  
}
  


#' Surface conductance to water vapor
#' 
#' @description Calculates surface conductance to water vapor from the inverted Penman-Monteith
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
#' @param PM        Calculate Gs based on the inverted Penman-Monteith (PM) equation? 
#'                  Defaults to TRUE. If FALSE, a simple gradient-conductance approach is used.
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  Mw - molar mass of water vapor (kg mol-1)
#' 
#' 
#' @details Surface conductance (Gs) is calculated from the inverted Penman-Monteith equation:
#' 
#'  \deqn{Gs = ( LE * Ga * \gamma ) / ( \Delta * (Rn-G-S) + \rho * cp * Ga * VPD - LE * ( \Delta + \gamma ) )}
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
#'  If PM is set to FALSE, Gs is calculated from VPD and ET only:
#'  
#'  \deqn{Gs = ET/pressure * VPD}
#'  
#'  Note that this approach assumes fully coupled conditions (Ga = inf).
#'  
#' 
#' @return a dataframe with the following columns: 
#'  \item{Gs_ms}{Surface conductance in m s-1}
#'  \item{Gs_mol}{Surface conductance in mol m-2 s-1}
#' 
#' @examples 
#' load("DE-Tha_2010.rda") # load data
#' 
#' # filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_2010 <- filter.data(DE_Tha_2010,quality.control=TRUE,quality.ext="_qc",
#'                            vars.qc=c("precip","Tair","VPD","GPP_nt","H","LE"),
#'                            tprecip=0.01,precip.hours=24,NA.as.precip=F,trad=200,
#'                            ttemp=5,tustar=0.2,tGPP=0.5,ws=15,min.int=5,trH=0.95) 
#' 
#' # select data that were filtered as in the previous function and in June 
#' DE_Tha_June_2010 <- DE_Tha_2010[DE_Tha_2010[,"valid"] == 1 & DE_Tha_2010[,"month"] == 6,]
#' 
#' # calculate Gs based on a simple gradient approach
#' Gs_gradient <- surface.conductance(DE_Tha_June_2010,Tair="Tair",pressure="pressure",VPD="VPD",PM=FALSE)
#' summary(Gs_gradient)
#' 
#' # calculate Gs from the the inverted PM equation (now Rn, and Ga are needed)
#' # calculate a simple estimate of Ga based on Thom 1972
#' Ga <- aerodynamic.conductance(DE_Tha_June_2010,stab_correction=F,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default). Note that
#' # note that Ga is not added to the data.frame 'DE_Tha_June_2010'
#' Gs_PM <- surface.conductance(DE_Tha_June_2010,Tair="Tair",pressure="pressure",Rn="Rn",G="G",S=NULL,
#'                              VPD="VPD",Ga=Ga,PM=TRUE)
#' summary(Gs_PM)
#' 
#' 
#' # compare to an Ga estimate that accounts for stability effects
#' Ga_stab <- aerodynamic.conductance(DE_Tha_June_2010,zr=42,zh=26.5,d=18.5,z0m=2.6,stab_correction=T,
#'                               Rb_model="Thom_1972")[,"Ga_h"]
#'                               
#' # now add Ga2 to the data.frame 'DE_Tha_June_2010'
#' DE_Tha_June_2010$Ga_stab <- Ga_stab
#' Gs_PM2 <- surface.conductance(DE_Tha_June_2010,Tair="Tair",pressure="pressure",Rn="Rn",G="G",S=NULL,
#'                              VPD="VPD",Ga="Ga_stab",PM=TRUE)
#' summary(Gs_PM2)
#'                              
#' @export
surface.conductance <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga",missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                PM=TRUE,constants=bigleaf.constants()){ 
    Tair     <- check.columns(data,Tair)
    pressure <- check.columns(data,pressure)
    VPD      <- check.columns(data,VPD)
    LE       <- check.columns(data,LE)
     
   if (!PM){
   
    Gs_simple_mol <- (LE.to.ET(LE,Tair)/constants$Mw) * pressure / VPD
    Gs_simple_ms  <- mol.to.ms(Gs_simple_mol,Tair,pressure)
     
    return(data.frame(Gs_ms=Gs_simple_ms,Gs_mol=Gs_simple_mol))
   
   } else {
  
    Rn       <- check.columns(data,Rn)
    Ga       <- check.columns(data,Ga)

      
     if(!is.null(G)){
       G <- check.columns(data,G)
       if (!missing.G.as.NA){G[is.na(G)] <- 0}
     } else {
       warning("Ground heat flux G is not provided and set to 0.")
       G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
     }
     
     if(!is.null(S)){
       S <- check.columns(data,S)
       if(!missing.S.as.NA){S[is.na(S)] <- 0 }
     } else {
       warning("Energy storage fluxes S are not provided and set to 0.")
       S <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
     }
     
    
     Delta <- Esat(Tair)[,"Delta"]
     gamma <- psychrometric.constant(Tair,pressure)
     rho   <- air.density(Tair,pressure)
   
     Gs_ms  <- ( LE * Ga * gamma ) / ( Delta * (Rn-G-S) + rho * constants$cp * Ga * VPD - LE * ( Delta + gamma ) )
     Gs_mol <- ms.to.mol(Gs_ms,Tair,pressure)
     
     return(data.frame(Gs_ms,Gs_mol))
   
   } 

}






#' Big-leaf surface conditions
#' 
#' @description Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' @param data       Data.frame or matrix containing all required input variables
#' @param Tair       Air temperature (deg C)
#' @param pressure   Atmospheric pressure (kPa)
#' @param Rn         Net radiation (W m-2)
#' @param H          Sensible heat flux (W m-2)
#' @param LE         Latent heat flux (W m-2)
#' @param VPD        Vapor pressure deficit (kPa)
#' @param Ga         Aerodynamic conductance for heat and water vapor (m s-1)
#' @param calc.Csurf Calculate surface CO2 concentration?
#' @param Ca         Atmospheric CO2 concentration (mol mol-1). Required if calc.Ca = TRUE
#' @param NEE        Net ecosystem exchange (umol m-2 s-1). Required if calc.Ca = TRUE
#' @param Ga_CO2     Aerodynamic conductance for CO2 (mol m-2 s-1). Required if calc.Ca = TRUE          
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
#' @examples 
#' # load("FR-Pue_2005.rda") # load data for the site FR-Pue  
#' # select data from July
#' FR_Pue_2005_July <- FR_Pue_2005[FR_Pue_2005[,"month"] == 7,]
#'   
#'       
#' # calculate a simple estimate of Ga based on Thom 1972
#' Ga <- aerodynamic.conductance(FR_Pue_2005_July,stab_correction=F,Rb_model="Thom_1972")[,"Ga_h"]
#'      
#' # calculate surface conditions (without Ca)      
#' surf <- surface.conditions(FR_Pue_2005_July,Ga=Ga,calc.Ca=F)           
#' summary(surf)  # note the outliers in both directions. They are associated with very low values of Ga    
#' 
#' # calculate surface conditions including Ca at the surface (NEE and Ga for CO2 are now needed as input)                                                                            
#' # Ga for CO2 can be calculated with the same function as above:
#' Ga_CO2_ms <- aerodynamic.conductance(FR_Pue_2005_July,stab_correction=F,Rb_model="Thom_1972")[,"Ga_CO2"]
#' # note that the function requires Ga for CO2 in mol m-2 s-1, not in ms-1.
#' Ga_CO2_mol <- ms.to.mol(Ga_CO2_ms,FR_Pue_2005_July[,"Tair"],FR_Pue_2005_July[,"pressure"])
#' surf <- surface.conditions(FR_Pue_2005_July,Ga=Ga_CO2_mol,NEE="NEE",calc.Ca=T)                                                                                        
#'                                                                                                                                                                                                                                                                                                            
#' @export 
surface.conditions <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",Ca="Ca",
                               VPD="VPD",Ga="Ga",calc.Csurf=F,Ga_CO2="Ga_CO2",NEE="NEE",
                               constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  LE       <- check.columns(data,LE)
  H        <- check.columns(data,H)
  VPD      <- check.columns(data,VPD)
  Ga       <- check.columns(data,Ga)
  
  if (calc.Csurf){
    Ca       <- check.columns(data,Ca)
    NEE      <- check.columns(data,NEE)
    Ga_CO2   <- check.columns(data,Ga_CO2)
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
    Ca_surf <- Ca.surface(Ca,NEE,Ga_CO2)
  }
  
  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
}



#' Canopy surface CO2 concentration
#'
#' @description the CO2 concentration at the canopy surface derived from net ecosystem
#'              CO2 exchange and measured atmospheric CO2 concentration.
#'              
#' @param Ca       atmospheric CO2 concentration (umol mol-1)
#' @param NEE      net ecosystem exchange (umol CO2 m-2 s-1)
#' @param Ga_CO2   aerodynamic conductance for CO2 (mol m-2 s-1)
#' 
#' @details CO2 concentration at the canopy surface is calculated as:
#' 
#'        \deqn{Ca_surf = Ca + NEE / Ga_CO2}
#'        
#'        Note that this equation can be used for any gas measured (with NEE
#'        replaced by the net exchange of the gas and Ga_CO2 by the Ga of 
#'        the respective gas).
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
#' # determine Trad of an object that has an emissivity of 0.98 and emits longwave radiation of 400Wm-2  
#' Trad.surface(400,0.98)
#' 
#' @export
Trad.surface<- function(longwave.up,emissivity,constants=bigleaf.constants()){
  
  Trad.K    <- (longwave.up / (constants$sigma * emissivity))^(1/4)
  Trad.degC <- Trad.K - constants$Kelvin
  
  return(data.frame(Trad.K,Trad.degC))
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
#'          canopy):
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1}{\epsilon + 1 + \frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 1) / ( \epsilon + 1 + Ga/Gc)}
#'          
#'          where \eqn{\epsilon = \frac{s}{\gamma}}{\epsilon = s/\gamma} is a dimensionless coefficient
#'          with s being the slope of the saturation vapor pressure curve (Pa K-1), and \eqn{\gamma} the 
#'          psychrometric constant (Pa K-1).
#'          
#'          The approach "Martin_1989" by Martin 1989 additionally takes radiative coupling
#'          into account:
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1 + \frac{Gr}{Ga}}{\epsilon + (1 + \frac{Ga}{Gs}) (1 + \frac{Gr}{Ga})}}{%
#'          \Omega = (\epsilon + 1 + Gr/Ga) / (\epsilon + (1 + Ga/Gs) (1 + Gr/Ga))}
#' 
#' @return omega The decoupling coefficient omega (-)
#' 
#' @references Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49. 
#'             
#'             Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#' 
#' @export
decoupling <- function(data,Tair="Tair",pressure="pressure",Ga="Ga",Gs="Gs",
                       approach=c("JarvisMcNaughton_1986","Martin_1989"),
                       LAI,constants=bigleaf.constants()){

  approach    <- match.arg(approach)
  canopy.type <- match.arg(canopy.type)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  Ga       <- check.columns(data,Ga)
  Gs       <- check.columns(data,Gs)

  Delta   <- Esat(Tair)[,"Delta"]
  gamma   <- psychrometric.constant(Tair,pressure,constants=constants)
  epsilon <- Delta/gamma
  
  if (approach == "JarvisMcNaughton_1986"){
    
    omega <- (epsilon + 1) / (epsilon + 1 + Ga/Gs)
    
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



#' Wet bulb temperature
#' 
#' @description calculates the wet bulb temperature, i.e. the temperature
#'              that the air would have if it was saturated.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' @details Wet-bulb temperature (Tw) is calculated from the following expression:
#'          
#'          \deqn{e = Esat(Tw) - gamma* (Tair - Tw)}
#'          
#'          The equation is solved for Tw using \link{\code{nls}}.
#'          Actual vapor pressure e is calculated from VPD using the function \code{\link{VPD.to.e}}.
#'          
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'             
#' @examples 
#' wetbulb.temp(Tair=c(20,25),pressure=100,VPD=c(1,1.6))             
#'             
#' @export
wetbulb.temp <- function(Tair,pressure,VPD,constants=bigleaf.constants()){

  gamma  <- psychrometric.constant(Tair,pressure)
  ea     <- VPD.to.e(VPD,pressure,Tair)
  
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
#' @param Tair      Air temperature (degC)
#' @param presssure Air pressure (kPa)                           
#' @param VPD       Vapor pressure deficit (kPa)
#' 
#' @details Dew point temperature (Td) is defined by:
#' 
#'         \deqn{e = Esat(Td)}
#'    
#'          which is solved for Td using \link{\code{nls}}.
#'          
#' @return dew point temperature Td (deg C)
#' 
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#' 
#' @export              
dew.point <- function(Tair,pressure,VPD){
  
  ea     <- VPD.to.e(VPD,pressure,Tair)
  
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
#' @param pressure  Air pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' @details the virtual temperature is given by:
#'  
#'  \deqn{Tv = Tair / (1 - (1 - eps) e/pressure)}
#' 
#'  where Tair is in Kelvin (converted internally). Likewise, VPD is converted 
#'  to actual vapor pressure (e in kPa) with \link{\code{VPD.to.e}} internally.
#' 
#' @return virtual temperature (deg C)
#' 
#' @references Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'  
#' @examples 
#' virtual.temp(25,100,1.5)                        
#'               
#' @export
virtual.temp <- function(Tair,pressure,VPD,constants=bigleaf.constants()){
   
  e    <- VPD.to.e(VPD,pressure,Tair)
   Tair <- Tair + constants$Kelvin
   
   Tv <- Tair / (1 - (1 - constants$eps) * e/pressure) 
   Tv <- Tv - constants$Kelvin
   
   return(Tv)
}




#' Saturation vapor pressure (Esat) and slope of the Esat curve
#' 
#' @description Calculates Saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' @param Tair     Air temperature (deg C)
#' @param formula  Formula to be used. Either "Sonntag_1990" (Default) or "Alduchov_1996"
#' 
#' @details Esat (kPa) is calculated based on the Magnus equation:
#' 
#'  \deqn{Esat = a * exp((b * Tair) / (c + Tair)) / 1000}
#'  
#'  where the coefficients a, b, c take different values depending on the formula used.
#'  The slope of the Esat curve (\eqn{\Delta}) is calculated as the first derivative of the function:
#'  
#'  formula
#' 
#' @return A dataframe with the following columns: 
#'    \item{Esat}{Saturation vapor pressure (kPa)}
#'    \item{Delta}{Slope of the saturation vapor pressure curve (kPa K-1)}
#'    
#'    
#' @references Sonntag_1990 Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90, and psychrometric formulae 
#'             
#'             Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure. Journal of Applied Meteorology, 1996, 35, 601-609
#' 
#' @examples 
#' Esat(seq(0,45,5))[,"Esat"]  # Esat in kPa
#' Esat(seq(0,45,5))[,"Delta"] # the corresponding slope of the Esat curve in kPa K-1
#'             
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





#' Conversions between global radiation and photosynthetic photon flux density
#' 
#' @description Converts radiation from W m-2 to umol m-2 s-1 and vice versa.
#' 
#' @param Rg       Global radiation = incoming short-wave radiation at the surface (W m-2)
#' @param PPFD     Photosynthetic photon flux density (umol m-2 s-1)
#' @param J_to_mol Conversion factor from J m-2 s-1 (=W m-2) to umol (quanta) m-2 s-1
#' @param frac_PAR Fraction of incoming solar irradiance that is photosynthetically 
#'                 active radiation (PAR); defaults to 0.5
#'                 
#' @details 
#' The conversion is given by:
#' 
#'  \deqn{PPFD = Rg * frac_PAR * J_to_mol}
#'  
#' by default, the combined conversion factor is 2.3
#'
#' @examples 
#' # convert a measured incoming short-wave radiation of 500 Wm-2 to 
#' # PPFD in umol m-2 s-1 and backwards
#' Rg.to.PPFD(500)
#' PPFD.to.Rg(1150)
#' 
#' @export
Rg.to.PPFD <- function(Rg,J_to_mol=4.6,frac_PAR=0.5){
  PPFD <- Rg * frac_PAR * J_to_mol
  return(PPFD)
}


#' @rdname Rg.to.PPFD
#' @export
PPFD.to.Rg <- function(PPFD,J_to_mol=4.6,frac_PAR=0.5){
  Rg <- PPFD / frac_PAR / J_to_mol
  return(Rg)
} 




#' Conversion between mass and molar units
#' 
#' @description Converts CO2 quantities from umol CO2 m-2 s-1 to gC and vice versa.
#' 
#' @param CO2_flux  CO2 flux (umol CO2 m-2 s-1)
#' @param C_flux    C flux (gC m-2 s-1)
#' @param constants Cmol - molar mass of carbon (kg mol-1)
#' 
#' @details The conversions are given by:
#' 
#'         
#' @export
umolCO2.to.gC <- function(CO2_flux,constants=bigleaf.constants()){
  
  C_flux <- CO2_flux * 1e-06 * constants$Cmol * 1000
  
  return(C_flux)
}


#' @rdname umolCO2.to.gC
#' @param C_flux CO2 flux (gC m-2 s-1)
#'
gC.to.umolCO2 <- function(C_flux,constants=bigleaf.constants()){

  CO2_flux <- C_flux / 1000 * (constants$Cmol / 2*constants$Omol) * 1e06

  return(CO2_flux)
}



###########################
#### Evapotranspiration ###-------------------------------------------------------------------------------
###########################

#' Potential evapotranspiration
#' 
#' @description Potential evapotranspiration ET_pot according to Priestley & Taylor 1972.
#' 
#' @param data      Data.frame or matrix containing all required variables; optional
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param alpha     Priestley-Taylor coefficient (-)
#' @param missing.G.as.NA  if TRUE, missing G are treated as NA,otherwise set to 0. 
#' @param missing.S.as.NA  if TRUE, missing S are treated as NA,otherwise set to 0. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @details Potential evapotranspiration is calculated according to Priestley & Taylor, 1972:
#' 
#'          \deqn{LE_pot = (\alpha * \Delta * (Rn - G - S)) / (\Delta + \gamma)}
#'
#' @return a data.frame with the following columns:
#'         \item{ET_pot}{potential evapotranspiration (kg m-2 s-1)}
#'         \item{LE_pot}{potential latent heat flux (W m-2)}
#'         
#' @note If the first argument 'data' is provided (either a matrix or a data.frame),
#'       the following variables can be provided as character (in which case they are interpreted as
#'       the column name of 'data') or as numeric vectors, in which case they are taken
#'       directly for the calculations. If 'data' is not provided, all input variables have to be
#'       numeric vectors.        
#'   
#' @references Priestley C.H.B., Taylor R.J., 1972: On the assessment of surface heat flux
#'             and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
#'             
#' @examples 
#' # Calculate potential ET from a surface that receives Rn of 400 Wm-2
#' ET.pot(Tair=30,pressure=100,Rn=400,alpha=1.26)    
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
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat(Tair)[,"Delta"]
  
  LE_pot <- (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
  ET_pot <- LE.to.ET(LE_pot,Tair)
  
  return(data.frame(ET_pot,LE_pot))
}




#' Reference evapotranspiration
#' 
#' @description Reference evapotranspiration calculated from the Penman-Monteith
#'              equation with a pre-defined surface conductance.
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param Gs        Surface conductance (m s-1); defaults to 0.0143 ms-1 (~ 0.58mol m-2 s-1)
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if TRUE, missing G are treated as NA,otherwise set to 0. 
#' @param missing.S.as.NA  if TRUE, missing S are treated as NA,otherwise set to 0. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                                
#' @details Reference evapotranspiration is calculated according to the Penman-Monteith
#'          equation:
#' 
#'          \deqn{LE_0 = (\Delta * (Rn - G - S) * \rho * cp * VPD * Ga) / (\Delta + \gamma * (1 + Ga/Gs)}
#'          
#'          The reference evapotranspiration is calculated with respect to a 'reference surface',
#'          which is typically a well-watered grass/crop of 0.12 height, an albedo of 0.23 and 
#'          a surface conductance of ~ 0.6 mol m-2 s-1 (Allen et al. 1998), but can be calculated for any other
#'          surface.
#'
#' @return a data.frame with the following columns:
#'         \item{ET_0}{reference evapotranspiration (kg m-2 s-1)}
#'         \item{LE_0}{reference latent heat flux (W m-2)}              
#'                  
#' @references  Allen R.G., Pereira L.S., Raes D., Smith M., 1998: Crop evapotranspiration -
#'              Guidelines for computing crop water requirements - FAO Irrigation and drainage
#'              paper 56.
#' 
#' @examples 
#' # Calculate ET_ref for a surface with known Gs (=0.5mol m-2 s-1) and Ga (0.1 ms-1)
#' 
#' # first, Gs has to be converted to ms-1
#' Gs_ms <- mol.to.ms(0.5,Tair=20,pressure=100)
#' ET_ref <- ET.ref(Gs=Gs_ms,Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400)
#' 
#' # now cross-check with the inverted version
#' surface.conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=ET_ref[,"LE_ref"])
#' 
#' @export                 
ET.ref <- function(data,Gs=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga",
                   G=NULL,S=NULL,missing.G.as.NA=F,missing.S.as.NA=F,
                   constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  VPD      <- check.columns(data,VPD)
  Rn       <- check.columns(data,Rn)
  Ga       <- check.columns(data,Ga)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat(Tair)[,"Delta"]
  rho    <- air.density(Tair,pressure)
  
  LE_ref <- (Delta * (Rn - G - S) + rho * constants$cp * VPD * Ga) / 
          (Delta + gamma * (1 + Ga / Gs))
  
  ET_ref <- LE.to.ET(LE_ref,Tair)
  
  return(data.frame(ET_ref,LE_ref))
  
}





#' Imposed and Equilibrium Evapotranspiration
#' 
#' @description Evapotranspiration (ET) split up into imposed ET and equilibrium ET.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param VPD       Air vapor pressure deficit (kPa)
#' @param Gs        surface conductance (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if TRUE, missing G are treated as NA,otherwise set to 0. 
#' @param missing.S.as.NA  if TRUE, missing S are treated as NA,otherwise set to 0.
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' @details Total evapotranspiration can be written in the form:
#' 
#'          \deqn{ET = \Omega ET_eq + (1 - \Omega)ET_imp}
#'          
#'          where ET_eq is the equilibrium evapotranspiration rate, the ET rate 
#'          that would occur under uncoupled conditions, where the heat budget
#'          is dominated by radiation (when Ga -> 0).
#'          
#'          \deqn{ET_eq = (\Delta * (Rn - G - S) * \lambda) / (\Delta + \gamma)}
#'          
#'          and ET_imp is the imposed evapotranspiration rate, the ET rate
#'          that would occur under fully coupled conditions (when Ga -> inf):
#'          
#'          \deqn{ET_imp = (\rho * cp * VPD * Gs * \lambda) / \gamma}
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
#' @export
ET.components <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs",
                          Rn="Rn",G=NULL,S=NULL,missing.G.as.NA=F,missing.S.as.NA=F,
                          constants=bigleaf.constants()){
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  VPD      <- check.columns(data,VPD)
  Rn       <- check.columns(data,Rn)
  Gs       <- check.columns(data,Gs)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  rho    <- air.density(Tair,pressure)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat(Tair)[,"Delta"]
  
  LE_eq  <- (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp <- (rho * constants$cp * Gs * VPD) / gamma
  
  ET_imp <- LE.to.ET(LE_imp,Tair)
  ET_eq  <- LE.to.ET(LE_eq,Tair)
  
  return(data.frame(ET_eq,ET_imp,LE_eq,LE_imp))
}





############################
#### Big-Leaf Physiology ###-------------------------------------------------------------------------------
############################

#' Bulk intercellular CO2 concentration
#' 
#' @description Bulk canopy intercellular CO2 concentration (Ci) calculated based on Fick's law
#'              given surface conductance (Gs), gross primary productivity (GPP) and 
#'              atmospheric CO2 concentration (Ca).
#'                            
#' @param Ca         Atmospheric CO2 concentration (umol mol-1)              
#' @param GPP        Gross primary productivity (umol CO2 m-2 s-1)
#' @param RecoLeaf   Ecosytem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0          
#' @param Gs         Surface conductance to water vapor (mol m-2 s-1)
#' @param calc.Csurf Should the derived surface CO2 concentration be used instead of 
#'                   measured atmospheric CO2? If TRUE, Ca is derived as shown in \code{Details}.
#' @param Ga         Aerodynamic conductance to CO2 (mol m-2 s-1) 
#' @param NEE        Net ecosystem exchange (umol CO2 m-2 s-1), negative values indicate CO2 uptake by the ecosytem
#' @param constants  DwDc - Ratio of the molecular diffusivities for water vapor and CO2 (-)
#' 
#' @details Bulk intercellular CO2 concentration (Ci) is given by:
#' 
#'          \deqn{Ci = Ca - (GPP - Reco_Leaf)/(Gs/1.6)}
#'          
#'          where Gs/1.6 (mol m-2 s-1) represents the surface conductance to CO2.
#'          Note that Gs is required in mol m-2 s-1 for water vapor. Gs is converted to
#'          its value for CO2 internally.
#'          
#' @note The equation is based on Fick's law of diffusion and is equivalent to the
#'       often used equation at leaf level (ci = ca - An/gs).
#'       Note that GPP and Gs have a different interpretation than An and gs.
#'       Gs comprises non-pyhsiological contributions (i.e. physical evaporation)
#'       and is confounded by physical factors (e.g. energy balance non-closure).
#'       GPP does not account for dark respiration and is further subject to uncertainties
#'       in the NEE partitioning algorithm used. Leaf respiration can be provided,
#'       but it is usually not known at ecosytem level (as a consequence, Ci is likely to be 
#'       slightly underestimated)
#'       This function should be used with care and the resulting Ci might not be
#'       readily comparable to its leaf-level analogue and/or physiological meaningful.          
#' 
#' @return Ci - Bulk canopy intercellular CO2 concentration (umol mol-1)
#' 
#' @examples 
#' # calculate bulk canopy Ci of a productive ecosystem
#' intercellular.CO2(Ca=400,GPP=40,Gs=0.7)
#' 
#' # now calculate bulk canopy Ci, but with Ca at the canopy surface (Ga and NEE are needed)
#' # Ga has to be in mol m-2 s-1 and for CO2, assuming that Ga for CO2 in ms-1 is 0.05
#' # alternatively, the function aerodynamic.conductance returns Ga for CO2
#' Ga_mol <- ms.to.mol(0.05,Tair=20,pressure=100)
#' 
#' intercellular.CO2(Ca=400,GPP=40,Gs=0.7,calc.Csurf=T,Ga=Ga_mol,NEE=-55) # note the sign convention for NEE
#' 
#' @references Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385.
#'             
#'             Keenan T., Sabate S., Gracia C., 2010: The importance of mesophyll conductance in
#'             regulating forest ecosystem productivity during drought periods. Global Change Biology
#'             16, 1019-1034.
#' 
#' @export
intercellular.CO2 <- function(Ca,GPP,RecoLeaf=NULL,Gs,calc.Csurf=F,Ga=NULL,NEE=NULL,
                              constants=bigleaf.constants()){
  
  
  if (calc.Csurf){
    Ca <- Ca.surface(Ca,NEE,Ga)
  }
  
  if (is.null(RecoLeaf)){
    RecoLeaf <- 0
    warning("Respiration from the leaves is ignored and set to 0.")
  }
  
  Ci <- Ca - (GPP - RecoLeaf)/(Gs/constants$DwDc)

  return(Ci)
  
}



#' Bulk canopy maximum carboxylation efficiency (Vcmax)
#' 
#' @description Bulk canopy maximum carboxylation rate (Vcmax) from bulk intercellular 
#'              CO2 concentration using the Farquhar et al. 1980 model for C3 photosynthesis.
#'              
#' @param Temp      Surface (or air) temperature (degC) 
#' @param GPP       Gross primary productivty (umol m-2 s-1)
#' @param Ci        Bulk canopy intercellular CO2 concentration (umol mol-1)
#' @param PPFD      Photosynthetic photon flux density (umol m-2 s-1) 
#' @param PPFD_sat  PPFD threshold, above which the canopy is considered to 
#'                  be light saturated.
#' @param Oi        intercellular O2 concentration (mol mol-1)
#' @param Kc25      Michaelis-Menten constant for CO2 at 25 degC (umol mol-1)
#' @param Ko25      Michaelis-Menten constant for O2 at 25 degC (mmol mol-1)
#' @param Kc_Ha     activation energy for Kc (kJ mol-1)
#' @param Ko_Ha     activation energy for Ko (kJ mol-1)
#' @param Vcmax_Ha  activation energy for Vcmax (kJ mol-1)
#' @param Vcmax_Hd  deactivation energy for Vcmax (kJ mol-1)
#' @param dS        entropy term (J mol-1 K-1)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  Rgas - universal gas constant (J mol-1 K-1)
#'                  
#' @details The maximum carboxylation rate (Vcmax) and Vcmax25 (Vcmax at 25degC) are calculated
#'          as at leaf level. The required variables Gs and Ci can be calculated from surface.conductance()
#'          and intercellular.CO2(), respectively.
#'          
#' @return a data.frame with the following columns:
#'         \item{Vcmax}{maximum carboxylation rate at Temp (umol m-2 s-1)}
#'         \item{Vcmax25}{maximum carboxylation rate at 25degC (umol m-2 s-1)}
#'         \item{Ko}{Michaelis-Menten constant for O2 at Temp (mmol mol-1)}
#'         \item{Kc}{Michaelis-Menten constant for CO2 at Temp (umol mol-1)}   
#'        
#' @references Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385. 
#'             
#'             Sharkey T.D., Bernacchi C.J., Farquhar G.D., Singsaas E.L., 2007: Fitting
#'             photosynthetic carbon dioxide response curves for C3 leaves. Plant, Cell and
#'             Environment 30, 1035-1040.
#'             
#'                             
#' @export                  
carbox.rate <- function(Temp,GPP,Ci,PPFD,PPFD_sat,Oi=0.21,Kc25=404.9,Ko25=278.4,
                        Kc_Ha=79.43,Ko_Ha=36.38,Vcmax_Ha=65.33,Vcmax_Hd=200,
                        dS=640,constants=bigleaf.constants()){
  
  Temp <- Temp + constants$Kelvin
  Tref <- constants$Kelvin + 25
  
  Kc_Ha    <- Kc_Ha * 1000
  Ko_Ha    <- Ko_Ha * 1000
  Vcmax_Ha <- Vcmax_Ha * 1000
  Vcmax_Hd <- Vcmax_Hd * 1000
  
  # Temperature dependencies of photosynthetic parameters 
  Kc <- Kc25 * exp(Kc_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  Ko <- Ko25 * exp(Ko_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  Ko <- Ko / 1000
  
  # exclude light unsaturated GPP 
  GPP[PPFD < PPFD_sat] <- NA
  GPPc <- GPP
  
  # calculate Vcmax and Vcmax25
  Vcmax   <- (GPPc * (Ci + Kc*(1.0 + Oi/Ko))) / Ci
  
  Vcmax25 <- Vcmax / 
             ( exp(Vcmax_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp)) *
               (1 + exp((Tref*dS - Vcmax_Hd) / (Tref * constants$Rgas))) /
               (1 + exp((Temp*dS - Vcmax_Hd) / (Temp * constants$Rgas)))
             )
  
  
  return(data.frame(Vcmax,Vcmax25,Ko=Ko*1000,Kc))
}



#' Stomatal slope parameter "g1"
#' 
#' @description Estimation of the intrinsic WUE metric "g1" (stomatal slope) 
#'              from non-linear regression.
#' 
#' @param data      Data.frame or matrix containing all required columns
#' @param Tair      Air temperature (deg C)
#' @param pressure  Air pressure (kPa)
#' @param GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' @param Gs        Surface conductance to water vapor (mol m-2 s-1)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ca        Atmospheric CO2 concentration (umol mol-1)
#' @param model     Stomatal model used. One of c("USO","Ball&Berry","Leuning")
#' @param nmin      Minimum number of data required to perform the fit; defaults to 40.
#' @param fitg0     Should g0 and g1 be fitted simultaneously? 
#' @param g0        Minimum stomatal conductance (mol m-2 s-1); ignored if \code{fitg0} is TRUE.
#' @param fitD0     Should D0 be fitted along with g1 (and g0 if fitg0 = TRUE)?; only used if \code{model} is "Leuning"
#' @param D0        Stomatal sensitivity parameter to VPD; only used if \code{model} is "Leuning" and fitD0 is FALSE
#' @param Gamma     Canopy CO2 compensation point (umol mol-1); only used if \code{model} is "Leuning", defaults to 50 umol mol-1.
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  Rgas - universal gas constant (J mol-1 K-1)
#' 
#' @details All stomatal models were developed at leaf-level, but its parameters 
#'          can also be estimated at ecosystem level (but be aware of caveats).
#'          
#'          The unified stomatal optimization (USO) model is given as (Medlyn et al. 2011):
#'      
#'          \deqn{gs = g0 + 1.6*(1.0 + g1/sqrt(VPD)) * GPP/Ca}
#'          
#'          The empirical model by Ball et al. 1987 is defined as:
#'          
#'          \deqn{gs = g0 + g1* ((An * rH) / Ca)}
#'          
#'          Leuning 1995 suggested a revised version of the Ball&Berry model:
#'          
#'          \deqn{gs = g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0))}
#'          
#'          The parameters in the models are estimated using non-linear regression (\code{nls()}).
#'          Alternatively to measured VPD and Ca (i.e. conditions at instrument height), conditions at 
#'          the big leaf surface can be provided. They can be calculated using \code{\link{surface.conditions}}.
#'          
#' 
#' @return an nls model object, containing information on the fitted parameters, their uncertainty range,
#'         model fit, etc.
#' 
#' @references Medlyn, B.E., et al., 2011: Reconciling the optimal and empirical approaches to
#'             modelling stomatal conductance. Global Change Biology 17, 2134-2144.
#'             
#'             Ball, T.J., Woodrow, I.E., Berry, J.A. 1987: A model predicting stomatal conductance
#'             and its contribution to the control of photosynthesis under different environmental conditions.
#'             In: Progress in Photosynthesis Research, edited by J.Biggins, pp. 221-224, Martinus Nijhoff Publishers,
#'             Dordrecht, Netherlands.
#'             
#'             Leuning, R., 1995: A critical appraisal of a combined stomatal-photosynthesis
#'             model for C3 plants. Plant, Cell and Environment 18, 339-355.
#' 
#' @examples 
#' # calculate g1 for the site DE-Tha using data from June 2010
#' load("DE_Tha_2010") # load data
#' 
#' # filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_2010 <- filter.data(DE_Tha_2010,quality.control=TRUE,quality.ext="_qc",
#'                            vars.qc=c("precip","Tair","VPD","GPP_nt","H","LE"),
#'                            tprecip=0.01,precip.hours=24,NA.as.precip=F,trad=200,
#'                            ttemp=5,tustar=0.2,tGPP=0.5,ws=15,min.int=5,trH=0.95) 
#' 
#' # select data only that were filtered in the previous function and in June
#' DE_Tha_June_2010 <- DE_Tha_2010[DE_Tha_2010[,"valid"] == 1 & DE_Tha_2010[,"month"] == 6,]
#' 
#' # calculate Gs from the the inverted PM equation (now Rn, and Ga are needed)
#' # calculate a simple estimate of Ga based on Thom 1972
#' Ga <- aerodynamic.conductance(DE_Tha_June_2010,stab_correction=F,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default). Note that
#' # note that Ga is not added to the data.frame 'DE_Tha_June_2010'
#' Gs_PM <- surface.conductance(DE_Tha_June_2010,Tair="Tair",pressure="pressure",Rn="Rn",G="G",S=NULL,
#'                              VPD="VPD",Ga=Ga,PM=TRUE)[,"Gs_mol"]
#'                              
#' ### Estimate the stomatal slope parameter g1 using the USO model
#' mod_USO <- stomatal.slope(DE_Tha_June_2010,model="USO",GPP="GPP_nt",Gs=Gs_PM,nmin=40,fitg0=FALSE)
#' 
#' ### Estimate the same parameter from the Ball&Berry model and prescribe g0
#' mod_BB <- stomatal.slope(DE_Tha_June_2010,model="Ball&Berry",GPP="GPP_nt",Gs=Gs_PM,g0=0.01,nmin=40,fitg0=FALSE)                         
#' 
#' ## same for the Leuning model, but this time estimate both g1 and g0 (but fix D0)
#' mod_Leu <- stomatal.slope(DE_Tha_June_2010,model="Leuning",GPP="GPP_nt",Gs=Gs_PM,nmin=40,fitg0=TRUE,D0=1.5,fitD0=FALSE)                         
#' 
#' @export 
stomatal.slope <- function(data,Tair="Tair",pressure="pressure",GPP="GPP_nt",Gs="Gs",VPD="VPD",Ca="Ca",
                           model=c("USO","Ball&Berry","Leuning"),nmin=40,fitg0=FALSE,g0=0,fitD0=FALSE,
                           D0=1.5,Gamma=50,constants=bigleaf.constants()){
  
  model <- match.arg(model)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  GPP      <- check.columns(data,GPP)
  Gs       <- check.columns(data,Gs)
  VPD      <- check.columns(data,VPD)
  Ca       <- check.columns(data,Ca)
  
  nr_data <- sum(!is.na(GPP) & !is.na(Gs) & !is.na(VPD) & !is.na(Ca))
  
  if (nr_data < nmin){
    stop("number of data is less than 'nmin'. g1 is not fitted to the data.")
  } else {
    
    if (model == "USO"){
      
      if (fitg0){
        mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3))
      } else {
        mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3))
      }
      
    } else if (model == "Leuning"){
      
      if (fitg0){
        if (fitD0){
          mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5))
        } else {
          mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9))
        }
      } else {
        if (fitD0){
          mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5))
        } else {
          mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9))
        }
      }
      
    } else if (model == "Ball&Berry"){
      
      rH <- VPD.to.rH(VPD,Tair)
      
      if (fitg0){
        mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9))
      } else {
        mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9))
      }
      
    }
    
  }
  
  return(mod)
  
}




########################
#### Energy balance ####-------------------------------------------------------------------------------
########################

#' Biochemical energy
#' 
#' @description Radiant energy absorbed in photosynthesis or heat release by respiration calculated
#'              from net ecosystem exchange of CO2 (NEE).  
#' 
#' @param NEE     Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param alpha   Energy taken up/released by photosynthesis/respiration (J umol-1)
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
#' @examples 
#' # Calculate biochemical energy taken up by the ecosystem with a measured NEE of -30umol CO2 m-2 s-1             
#' biochemical.energy(NEE=-30)            
#'            
#' @export 
biochemical.energy <- function(NEE,alpha=0.422){
  Sp <- -alpha*NEE
  return(Sp)
}




#' Energy balance closure
#' 
#' @description Calculates the degree of the energy balance non-closure for the entire timespan
#'              based on the ratio of two sums (energy balance ratio), and ordinary least squares (OLS).
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
#'         if \code{instantaneous} is TRUE, only \code{EBR} is returned
#' 
#' @references Wilson K., et al. 2002: Energy balance closure at FLUXNET sites.
#'             Agricultural and Forest Meteorology 113, 223-243.
#'
#' @examples 
#' load("DE-Tha_2010.rda") # load data
#' 
#' ## characterize energy balance closure for DE-Tha in 2010
#' energy.closure(DE_Tha_2010,instantaneous=F)
#' 
#' ## look at half-hourly closure 
#' energy.closure(DE_Tha_2010,instantaneous=T)
#' 
#' @export
energy.closure <- function(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",instantaneous=F,
                           missing.G.as.NA=F,missing.S.as.NA=F){
  
  Rn <- check.columns(data,Rn)
  LE <- check.columns(data,LE)
  H  <- check.columns(data,H)
  
  if(!is.null(G)){
    G <- check.columns(data,G)
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    S <- check.columns(data,S)
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.")
    S <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }

  if (!instantaneous){
    comp <- complete.cases(Rn,G,S,LE,H)
    n    <- sum(comp)
  
    EBR <- sum(LE[comp] + H[comp]) / sum(Rn[comp] - G[comp] - S[comp])
  
    emod <- lm((LE + H) ~ (Rn - G - S))
    intercept <- summary(emod)$coef[1,1]
    slope     <- summary(emod)$coef[2,1]
    r_squared <- summary(emod)$r.squared
    
    return(c("n"=n,"intercept"=round(intercept,3),"slope"=round(slope,3),"r^2"=round(r_squared,3),"EBR"=round(EBR,3)))
  
  } else {
    
    EBR <- (LE + H) /(Rn - G - S)
    
    return(EBR)
  }
}   




#' Isothermal net radiation
#' 
#' @description Calculates the isothermal net radiation, i.e. the net radiation 
#'              that the surface would receive if it had the same temperature than
#'              the air.
#'              
#' @param data       Data.frame or matrix containing all required variables
#' @param Rn         Net radiation (W m-2)
#' @param Tair       Air temperature (degC)
#' @param Tsurf      Surface temperature (degC)
#' @param emissivity Emissivity of the surface (-)
#' @param constants  sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                   Kelvin - conversion degree Celsius to Kelvin 
#'
#' @details The isothermal net radiation (Rni) is given by:
#'          
#'          \deqn{Rni = Rn + \epsilon * \sigma * (Tsurf^4 - Tair^4)}
#'          
#'          with Tsurf and Tair in Kelvin.
#'          
#' @return Rni isothermal net radiation (W m-2)
#' 
#' @references Jones, H. 2014: Plants and Microclimate. 3rd edition, Cambridge
#'             University Press.
#' 
#' @examples 
#' # calculate isothermal net radiation of a surface that is 2?c warmer than the air.
#' isothermal.Rn(Rn=400,Tair=25,Tsurf=27,emissivity=0.98)
#' 
#' @export
isothermal.Rn <- function(data,Rn="Rn",Tair="Tair",Tsurf="Tsurf",emissivity,
                          constants=bigleaf.constants()){
  
  if (!missing(data)){
    Rn    <- check.columns(data,Rn)
    Tair  <- check.columns(data,Tair)
    Tsurf <- check.columns(data,Tsurf)
  }
  
  Tair  <- Tair + constants$Kelvin
  Tsurf <- Tsurf + constants$Kelvin
  
  Rni <- Rn + emissivity * constants$sigma * (Tsurf^4 - Tair^4)
  

  return(Rni)
  
}
