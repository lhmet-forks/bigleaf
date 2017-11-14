########################
### Filter functions ###----------------------------------------------------------------------------
########################

#' Basic data filtering
#'
#' @description Filters timeseries of EC data for high-quality values and specified
#'              meteorological conditions.
#' 
#' @param data            Data.frame or matrix containing all required input variables in 
#'                        half-hourly or hourly resolution. Including year, month, day information
#' @param quality.control Should quality control be applied? Defaults to TRUE.
#' @param vars.qc         Character vector indicating the variables for which quality filter should 
#'                        be applied. Ignored if \code{quality.control} is FALSE.
#' @param filter.growseas Should data be filtered for growing season? Defaults to FALSE.
#' @param filter.precip   Should precipitation filtering be applied? Defaults to FALSE.
#' @param filter.vars     Additional variables to be filtered. Vector of type character.
#' @param filter.vals.min Minimum values of the variables to be filtered. Numeric vector of 
#'                        the same length than \code{filter.vars}. Set to NA to be ignored.
#' @param filter.vals.max Maximum values of the variables to be filtered. Numeric vector of 
#'                        the same length than \code{filter.vars}. Set to NA to be ignored.
#' @param NA.as.invalid   If TRUE (the default) missing data are filtered out (applies to all variables).
#' @param quality.ext     The extension to the variables' names that marks them as 
#'                        quality control variables. Ignored if \code{quality.control} is FALSE.                       
#' @param good.quality    Which values indicate good quality (i.e. not to be filtered) 
#'                        in the quality control (qc) variables? Ignored if \code{quality.control} is FALSE.
#' @param missing.qc.as.bad If quality control variable is NA, should the corresponding data point be
#'                          treated as bad quality? Defaults to TRUE. Ignored if \code{quality.control} is FALSE.                        
#' @param precip          Precipitation (mm time-1)
#' @param GPP             Gross primary productivity (umol m-2 s-1); Ignored if \code{filter.growseas} is FALSE.
#' @param doy             Day of year; Ignored if \code{filter.growseas} is FALSE.
#' @param year            Year; Ignored if \code{filter.growseas} is FALSE.
#' @param tGPP            GPP threshold (fraction of 95% quantile of the GPP time series).
#'                        Must be between 0 and 1. Ignored if \code{filter.growseas} is FALSE.
#' @param ws              Window size used for GPP time series smoothing. 
#'                        Ignored if \code{filter.growseas} is FALSE.
#' @param min.int         Minimum time interval in days for a given state of growing season.
#'                        Ignored if \code{filter.growseas} is FALSE.
#' @param tprecip         Precipitation threshold used to demark a precipitation event (mm). 
#'                        Ignored if \code{filter.precip} is FALSE
#' @param precip.hours    Number of hours removed following a precipitation event (h).
#'                        Ignored if \code{filter.precip} is FALSE
#' @param records.per.hour Number of observations per hour. I.e. 2 for half-hourly data.
#' 
#' 
#' @details This routine consists of two parts:
#' 
#'          1) Quality control: All variables included in \code{vars.qc} are filtered for 
#'             good quality data. For these variables a corresponding quality variable with 
#'             the same name as the variable plus the extension as specified in \code{quality.ext}
#'             must be provided. For timesteps where the value of the quality indicator is not included
#'             in the argument \code{good.quality}, i.e. the quality is not considered as 'good', 
#'             its value is set to NA.
#'             
#'          2) Meteorological filtering. Under certain conditions (e.g. low ustar), the assumptions
#'             of the EC method are not fulfilled. Further, some data analysis require certain meteorological
#'             conditions, such as periods without rainfall, or active vegetation (growing season, daytime).
#'             The filter applied in this second step serves to exclude time periods that do not fulfill the criteria
#'             specified in the arguments. More specifically, timeperiods where one of the variables is higher
#'             or lower than the specified thresholds (\code{filter.vals.min} and \code{filter.vals.max})
#'             are set to NA for all variables. If a threshold is set to NA, it will be ignored.
#'          
#' @return The same data frame as provided as first argument to the function
#'         with an additional column "valid", which indicates whether the data fulfill 
#'         the filtering criteria (1) or not (0).
#'         
#' @note Variables considered of bad quality (as specified by the corresponding quality control variables)      
#'       will be set to NA by this routine. Data that do not fulfill the filtering critera are not set to
#'       NA but the corresponding "valid" entry is set to 0.
#' 
#' @examples 
#' # Example of data filtering; data are for a month within the growing season,
#' # hence growing season is not filtered.
#' DE_Tha_Jun_2014_2 <- filter.data(DE_Tha_Jun_2014,quality.control=FALSE,
#'                                  vars.qc=c("Tair","precip","H","LE"),
#'                                  filter.growseas=FALSE,filter.precip=TRUE,
#'                                  filter.vars=c("Tair","PPFD","ustar"),
#'                                  filter.vals.min=c(5,200,0.2),
#'                                  filter.vals.max=c(NA,NA,NA),NA.as.invalid=TRUE,
#'                                  quality.ext="_qc",good.quality=c(0,1),
#'                                  missing.qc.as.bad=TRUE,GPP="GPP_nt",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
#'                                  tprecip=0.1,precip.hours=24,records.per.hour=2)
#'
#'  # note the additional column 'valid' in DE_Tha_Jun_2014_2.
#'  # To remove timesteps marked as filtered out (i.e. 0 values in column 'valid'):
#'  DE_Tha_Jun_2014_2[DE_Tha_Jun_2014_2["valid"] == 0,] <- NA
#'   
#'   
#' @importFrom stats aggregate
#' @export                     
filter.data <- function(data,quality.control=TRUE,vars.qc=NULL,filter.growseas=FALSE,
                        filter.precip=FALSE,filter.vars=NULL,
                        filter.vals.min,filter.vals.max,NA.as.invalid=TRUE,
                        quality.ext="_qc",good.quality=c(0,1),
                        missing.qc.as.bad=TRUE,GPP="GPP",doy="doy",
                        year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
                        tprecip=0.01,precip.hours=24,records.per.hour=2){
  
  
  ### I) Quality control filter
  if (quality.control){
    
    if (is.null(vars.qc)){
      stop("quality.control (qc) is TRUE, but no qc variables are provided!")
    }
    
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
      
      cat(var,": ",qc_invalid," data points (",qc_invalid_perc,"%) set to NA",fill=TRUE,sep="")
    }
    cat("-------------------------------------",fill=TRUE)
  }
  
  
  #### II) Data filter
  valid <- rep(1L,nrow(data))
  
  # 1) GPP
  growseas_invalid <- numeric()
  if(filter.growseas){
    doy    <- check.columns(data,doy)
    year   <- check.columns(data,year)
    date <- strptime(paste0(year,"-",doy),format="%Y-%j")
    GPP              <- check.columns(data,GPP)
    GPP_daily        <- aggregate(GPP,by=list(strftime(date)),mean,na.rm=TRUE)[,2]
    growing_season   <- filter.growing.season(GPP_daily,tGPP=tGPP,ws=ws,min.int=min.int)
    growseas_invalid <- which(sapply(growing_season,rep,48) == 0)
  }
  
  # 2) precipitation
  precip_invalid <- numeric()
  if (filter.precip){
    precip <- check.columns(data,precip)
    if (NA.as.invalid){
      precip_events <- which(precip > tprecip | is.na(precip))
    } else {
      precip_events <- which(precip > tprecip)
    }
    precip_invalid <- unique(as.numeric(unlist(sapply(precip_events, function(x) x:(min(x+precip.hours*records.per.hour,nrow(data),na.rm=TRUE))))))
  }
  
  # 3) all other filter variables (defined in filter.vars)
  invalids <- list(growseas_invalid,precip_invalid)
  
  if (!is.null(filter.vars)){
    for (var in filter.vars){
      v  <- which(filter.vars == var)
      vf <- v + 2
      assign(var,check.columns(data,var))
      if (NA.as.invalid){
        invalids[[vf]] <- which(get(var) < filter.vals.min[v] | get(var) > filter.vals.max[v] | is.na(get(var)))
      } else {
        invalids[[vf]] <- which(get(var) < filter.vals.min[v] | get(var) > filter.vals.max[v] & !is.na(get(var)))
      }
    }
  }
  
  # 4) calculate number and percentage of filtered values
  invalids_perc <- sapply(invalids, function(x) round((length(x)/nrow(data))*100,2))
  
  additional_invalids <- sapply(2:length(invalids), function(x) 
    length(setdiff(invalids[[x]],unique(unlist(invalids[1:(x-1)])))))
  
  additional_invalids_perc <- round(additional_invalids/nrow(data)*100,2)
  
  
  # 5) write to output
  var.names <- c("growing season","precipitation",filter.vars)
  cat("Data filtering:",fill=TRUE)
  
  cat(length(growseas_invalid)," data points (",invalids_perc[1],"%) excluded by growing season filter",fill=TRUE,sep="")
  
  invisible(sapply(c(1:(length(invalids)-1)), function(x) cat(additional_invalids[x]," additional data points (",
                                                              additional_invalids_perc[x],"%) excluded by ",var.names[x+1],
                                                              " filter (",length(unlist(invalids[x+1]))," data points = ",
                                                              invalids_perc[x+1]," % in total)",fill=TRUE,sep="")))
  
  
  invalid        <- unique(unlist(invalids))
  valid[invalid] <- 0
  
  excl_perc <- round((length(invalid)/nrow(data))*100,2)
  
  cat(length(invalid)," data points (",excl_perc,"%) excluded in total",fill=TRUE,sep="")
  cat(nrow(data) - length(invalid)," valid data points (",100-excl_perc,"%) remaining.",fill=TRUE,sep="")
  
  
  # 6) return input data frame with additional 'valid' column
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
#' @importFrom stats quantile filter                                 
#' @export  
filter.growing.season <- function(GPPd,tGPP,ws=15,min.int=5){
  
  if(sum(is.na(GPPd)) < 0.5*length(GPPd)){
    
    growseas      <- rep(1,length(GPPd))
    GPP_threshold <- quantile(GPPd,probs=0.95,na.rm=TRUE)*tGPP
    
    ## smooth GPP
    GPPd_smoothed <- filter(GPPd,method="convolution",filter=rep(1/ws,ws))
    
    ## set values at the beginning and end of the time series to the mean of the original values
    wsd <- floor(ws/2)
    GPPd_smoothed[1:wsd] <- mean(GPPd[1:(2*wsd)],na.rm=TRUE)
    GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] <- mean(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)],na.rm=TRUE)
    
    # check for occurence of missing values and set them to mean of the values surrounding them
    missing <- which(is.na(GPPd_smoothed))
    if (length(missing) > 0){
      if (length(missing) > 10){warning("Attention, there is a gap in 'GPPd' of length n = ",length(missing))}
      replace_val <- mean(GPPd_smoothed[max(1,missing[1] - 4):min((missing[length(missing)] + 4),length(GPPd_smoothed))],na.rm=TRUE)
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