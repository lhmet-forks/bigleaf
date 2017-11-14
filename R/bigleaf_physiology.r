############################
#### Big-Leaf Physiology ###-------------------------------------------------------------------------------
############################

#' Bulk intercellular CO2 concentration
#' 
#' @description Bulk canopy intercellular CO2 concentration (Ci) calculated based on Fick's law
#'              given surface conductance (Gs), gross primary productivity (GPP) and 
#'              atmospheric CO2 concentration (Ca).
#'                            
#' @param data       Data.Frame or matrix with all required columns                            
#' @param Ca         Atmospheric CO2 concentration (umol mol-1)              
#' @param GPP        Gross primary productivity (umol CO2 m-2 s-1)
#' @param Gs         Surface conductance to water vapor (mol m-2 s-1)
#' @param Reco_Leaf  Ecosytem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0          
#' @param calc.Csurf Should the derived surface CO2 concentration be used instead of 
#'                   measured atmospheric CO2? If TRUE, Ca is derived as shown in \code{Details}.
#' @param Ga_CO2     Aerodynamic conductance to CO2 (m s-1) 
#' @param NEE        Net ecosystem exchange (umol CO2 m-2 s-1), negative values indicate CO2 uptake by the ecosytem
#' @param Tair       Air temperature (degC); ignored if \code{calc.Csurf} is FALSE
#' @param pressure   Atmospheric pressure (kPa); ignored if \code{calc.Csurf} is FALSE
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
#' @return \item{Ci -}{Bulk canopy intercellular CO2 concentration (umol mol-1)}
#' 
#' @references Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385.
#'             
#'             Keenan T., Sabate S., Gracia C., 2010: The importance of mesophyll conductance in
#'             regulating forest ecosystem productivity during drought periods. Global Change Biology
#'             16, 1019-1034.
#'             
#' @examples 
#' # calculate bulk canopy Ci of a productive ecosystem
#' intercellular.CO2(Ca=400,GPP=40,Gs=0.7)
#' 
#' # now calculate bulk canopy Ci, but with Ca at the canopy surface (Ga and NEE are needed)
#' # The function \code{\link{aerodynamic.conductance}} can be used to calculate \code{Ga_CO2}.
#' # Here, \code{Ga_CO2} of 0.05ms-1 is assumed.
#' 
#' intercellular.CO2(Ca=400,GPP=40,Gs=0.7,calc.Csurf=TRUE,Ga_CO2=0.05,NEE=-55,
#'                   Tair=25,pressure=100) 
#' # note the sign convention for NEE
#' 
#' @export
intercellular.CO2 <- function(data,Ca,GPP,Gs,Reco_Leaf=NULL,calc.Csurf=FALSE,
                              Ga_CO2=NULL,NEE=NULL,Tair=NULL,pressure=NULL,
                              constants=bigleaf.constants()){
  
  check.input(data,list(Ca,GPP,Gs))

  if (calc.Csurf){
    
    check.input(data,list(Tair,pressure,Ga_CO2,NEE))
    
    Ca <- Ca.surface(Ca,NEE,Ga_CO2,Tair,pressure)
  }
  
  if (is.null(Reco_Leaf)){
    Reco_Leaf <- 0
    warning("Respiration from the leaves is ignored and set to 0.")
  }
  
  Ci <- Ca - (GPP - Reco_Leaf)/(Gs/constants$DwDc)
  
  return(Ci)
  
}



#' Bulk canopy biochemical parameters Vcmax and Jmax
#' 
#' @description Bulk canopy maximum carboxylation rate (Vcmax), electron transport rate (J),
#'              and maximum electron transport rate (Jmax) from bulk intercellular 
#'              CO2 concentration using the Farquhar et al. 1980 model for C3 photosynthesis.
#'           
#' @param data       Data.Frame or matrix with all required columns                 
#' @param Temp       Surface (or air) temperature (degC) 
#' @param GPP        Gross primary productivty (umol m-2 s-1)
#' @param Ci         Bulk canopy intercellular CO2 concentration (umol mol-1)
#' @param PPFD       Photosynthetic photon flux density (umol m-2 s-1) 
#' @param PPFD_j     PPFD threshold, below which the canopy is considered to 
#'                   be RuBP regeneration limited. Defaults to 500 umol m-2 s-1.
#' @param PPFD_c     PPFD threshold, above which the canopy is considered to 
#'                   be Rubisco limited. Defaults to 1000 umol m-2 s-1.
#' @param Reco_Leaf  Ecosytem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0 
#' @param Oi         Intercellular O2 concentration (mol mol-1)
#' @param Kc25       Michaelis-Menten constant for CO2 at 25 degC (umol mol-1)
#' @param Ko25       Michaelis-Menten constant for O2 at 25 degC (mmol mol-1)
#' @param Gam25      Photorespiratory CO2 compensation point ('Gamma star') 
#'                   at 25 degC (umol mol-1)
#' @param Kc_Ha      Activation energy for Kc (kJ mol-1)
#' @param Ko_Ha      Activation energy for Ko (kJ mol-1)
#' @param Gam_Ha     Activation energy for Gam (kJ mol-1)
#' @param Vcmax_Ha   Activation energy for Vcmax (kJ mol-1)
#' @param Vcmax_Hd   Deactivation energy for Vcmax (kJ mol-1)
#' @param Vcmax_dS   Entropy term for Vcmax (J mol-1 K-1)
#' @param Jmax_Ha    Activation energy for Jmax (kJ mol-1)
#' @param Jmax_Hd    Deactivation energy for Jmax (kJ mol-1)
#' @param Jmax_dS    Entropy term for Jmax (J mol-1 K-1)
#' @param Theta      Curvature term in the light response function of J (-)
#' @param alpha_canopy Canopy absorptance (-)
#' @param constants    Kelvin - conversion degree Celsius to Kelvin \cr
#'                     Rgas - universal gas constant (J mol-1 K-1)
#'                  
#' @details The maximum carboxylation rate (Vcmax) and Vcmax25 (Vcmax at 25degC) are calculated
#'          as at leaf level. The required variables Gs and Ci can be calculated from 
#'          \code{\link{surface.conductance}} and \code{\link{intercellular.CO2}}, respectively.
#'          
#'          Gas exchange parameters are taken from Bernacchi et al. 2001 (assuming
#'          an infinite mesophyll conductance)
#'          
#'          Vcmax is calculated from the photosynthesis model by Farquhar et al. 1980,
#'          assuming that net photosynthesis is Rubisco-limited (RuBP-satured carboxylation
#'          rate, i.e. light has to be saturating):
#'          
#'          \deqn{Vcmax = (GPP * (Ci + Kc*(1.0 + Oi/Ko))) / Ci}
#'          
#'          Vcmax at canopy level is assumed to follow the same temperature response
#'          as at leaf level. Hence, Vcmax at 25degC (Vcmax25) is calculated as 
#'          (see e.g. Kattge & Knorr 2007):
#'          
#'          \deqn{Vcmax25 = Vcmax / 
#'                         ( exp(Vcmax_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp)) *
#'                         (1 + exp((Tref*dS - Vcmax_Hd) / (Tref * constants$Rgas))) /
#'                         (1 + exp((Temp*dS - Vcmax_Hd) / (Temp * constants$Rgas)))
#'                         )
#'                }
#'          
#' @note   The critical assumption is that bulk canopy photosynthesis is Rubisco-
#'         limited. Therefore, only GPP data above a certain light threshold are
#'         used for the calculation of Vcmax. This threshold (\code{PPFD_c}) is
#'         likely dependent on the canopy of the ecosystem and presumably higher
#'         for dense than for open canopies. A threshold of 500 umol m-2 s-1 PPFD
#'         is a reasonable working assumption (see Kosugi et al. 2013).
#'         The photosynthesis parameters are also likely dependent on the ecosystem
#'         and growth conditions.         
#'          
#' @return a data.frame with the following columns:
#'         \item{Vcmax25}{maximum carboxylation rate at 25degC (umol m-2 s-1)}
#'         \item{Jmax25}{maximum electron transport rate at 25degC (umol m-2 s-1)}
#'        
#' @references Lloyd J. et al., 1995: A simple calibrated model of Amazon rainforest productivity
#'             based on leaf biochemical properties. Plant, Cell and Environment 18, 1129-1145.
#' 
#'             Rayment M.B., Loustau D., Jarvis P.G., 2002: Photosynthesis and respiration
#'             of black spruce at three organizational scales: shoot, branch and canopy.
#'             Tree Physiology 22, 219-229.
#' 
#'             Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385. 
#'             
#'             Ueyama M. et al, 2016: Optimization of a biochemical model with eddy covariance
#'             measurements in black spruce forests of Alaska for estimating CO2 fertilization
#'             effects. Agricultural and Forest Meteorology 222, 98-111.
#'             
#'             Bernacchi C.J., Singsaas E.L., Pimentel C., Portis JR A.R., Long S.P., 2001:
#'             Improved temperature response functions for models of Rubisco-limited
#'             photosynthesis. Plant, Cell and Environment 24, 253-259. 
#'             
#'             Bernacchi C.J., Pimentel C., Long S.P., 2003: In vivo temperature response
#'             functions of parameters required to model RuBP-limited photosynthesis.
#'             Plant, Cell and Environment 26, 1419-1430.
#'             
#'             Kattge J., Knorr W., 2007: Temperature acclimation in a biochemical
#'             model of photosynthesis: a reanalysis of data from 36 species.
#'             Plant, Cell and Environment 30, 1176-1190.
#'             
#'
#' @examples 
#' bigleaf.Vcmax.Jmax(Temp=20,GPP=15,Ci=300,PPFD=1000,PPFD_j=c(100,400),PPFD_c=1000)  
#'                                       
#' @export                  
bigleaf.Vcmax.Jmax <- function(data,Temp,GPP,Ci,PPFD,PPFD_j=c(100,400),PPFD_c=1000,
                               Reco_Leaf=NULL,Oi=0.21,Kc25=404.9,Ko25=278.4,Gam25=42.75,
                               Kc_Ha=79.43,Ko_Ha=36.38,Gam_Ha=37.83,Vcmax_Ha=65.33,Vcmax_Hd=200,
                               Vcmax_dS=635,Jmax_Ha=43.9,Jmax_Hd=200,Jmax_dS=640,
                               Theta=0.7,alpha_canopy=0.8,
                               constants=bigleaf.constants()){
  
  check.input(data,list(Temp,GPP,Ci,PPFD))
  
  Temp <- Temp + constants$Kelvin
  Tref <- 25.0 + constants$Kelvin
  
  Kc_Ha    <- Kc_Ha * 1000
  Ko_Ha    <- Ko_Ha * 1000
  Gam_Ha   <- Gam_Ha * 1000
  Vcmax_Ha <- Vcmax_Ha * 1000
  Vcmax_Hd <- Vcmax_Hd * 1000
  Jmax_Ha  <- Jmax_Ha * 1000
  Jmax_Hd  <- Jmax_Hd * 1000
  
  # Temperature dependencies of photosynthetic parameters 
  Kc  <- Kc25 * exp(Kc_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  Ko  <- Ko25 * exp(Ko_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  Gam <- Gam25 * exp(Gam_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  Ko <- Ko / 1000
  
  # Presumed limitation states 
  GPPc = GPPj <- GPP
  GPPj[PPFD < PPFD_j[1] | PPFD > PPFD_j[2] | is.na(PPFD)] <- NA
  GPPc[PPFD < PPFD_c | is.na(PPFD)] <- NA
  
  if (is.null(Reco_Leaf)){
    Reco_Leaf <- 0
    warning("Respiration from the leaves is ignored and set to 0.")
  }
  
  # calculate Vcmax and J (electron transport rate)
  Vcmax <- (GPPc-Reco_Leaf) * (Ci + Kc*(1.0 + Oi/Ko)) / (Ci - Gam)
  J     <- (GPPj-Reco_Leaf) * (4.0 * Ci + 8.0 * Gam) / (Ci - Gam)

  # calculate Jmax from J
  APPFD_PSII <- PPFD * alpha_canopy * 0.85 * 0.5
  dat <- data.frame(J,APPFD_PSII)
  dat <- dat[complete.cases(dat),]
  if (nrow(dat) > 0){
    Jmax <- sapply(c(1:nrow(dat)), function(x) tryCatch(summary(nls(J ~ c((APPFD_PSII + Jmax - 
                                                            sqrt((APPFD_PSII + Jmax)^2 - 
                                                            4.0 * Theta * APPFD_PSII * Jmax)) /
                                                            (2.0 * Theta)),
                                                            start=list(Jmax=50),data=dat[x,],
                                                            algorithm="port"))$coef[1],
                                                            error=function(err){NA}
                                                            )
                       )
  } else {
    warning("Not enough observations to calculate Jmax!")
    Jmax <- NA
  }
  
  # calculate Vcmax25 and Jmax25
  Vcmax25 <- Arrhenius.temp.response(Vcmax,Temp-constants$Kelvin,Ha=Vcmax_Ha/1000,
                                     Hd=Vcmax_Hd/1000,dS=Vcmax_dS,constants=constants)
  
  Jmax25 <- Arrhenius.temp.response(Jmax,Temp-constants$Kelvin,Ha=Jmax_Ha/1000,
                                    Hd=Jmax_Hd/1000,dS=Jmax_dS,constants=constants)
  
  Vcmax25_MEAN <- mean(Vcmax25,na.rm=T)
  Vcmax25_SE   <- sd(Vcmax25)/sqrt((sum(!is.na(Vcmax25))))
  Jmax25_MEAN  <- mean(Jmax25,na.rm=T)
  Jmax25_SE    <- sd(Jmax25)/sqrt((sum(!is.na(Jmax25))))
  
  return(c("Vcmax25_MEAN"=Vcmax25_MEAN,"Vcmax25_SE"=Vcmax25_SE,
           "Jmax25_MEAN"=Jmax25_MEAN,"Jmax25_SE"=Jmax25_SE))
}




#' (Modified) Arrhenius temperature response function
#' 
#' @description (Modified) Arrhenius function describing
#'              the temperature response of biochemical parameters.
#'              
#' @param param Parameter measured at measurement temperature (umol m-2 s-1)            
#' @param Temp  Measurement temperature (degC)
#' @param Ha    Activation energy for param (kJ mol-1)
#' @param Hd    Deactivation energy for param (kJ mol-1)
#' @param dS    Entropy term for param (J mol-1 K-1)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  Rgas - universal gas constant (J mol-1 K-1)
#'                  
#' @details The temperature response is given by a modified form of the Arrhenius
#'          function:
#' 
#'             \deqn{param25 = param / 
#'                         ( exp(Ha * (Temp - Tref) / (Tref*Rgas*Temp)) *
#'                         (1 + exp((Tref*dS - Hd) / (Tref * Rgas))) /
#'                         (1 + exp((Temp*dS - Hd) / (Temp * Rgas)))
#'                         )
#'                  }
#'                  
#'          where param is the value/rate of the parameter at measurement temperature,
#'          Temp is temperature in K, Tref is reference temperature (298.15K), and Rgas
#'          is the universal gas constant (8.314 J K-1 mol-1). Ha is the activation
#'          energy (kJ mol-1), Hd is the deactivation energy (kJ mol-1), and dS the
#'          entropy term (J mol-1 K-1) of the respective parameter.
#'          
#'          If either Hd or dS or both are not provided, the equation above reduces
#'          to the first term (i.e. the common Arrhenius equation without the deactivation
#'          term.)         
#'                                  
#' @return param25 - value of the input parameter at the reference temperature of 25degC (umol m-2 s-1)
#'               
#' @references Johnson F.H., Eyring H., Williams R.W. 1942: 
#'             The nature of enzyme inhibitions in bacterial luminescence: sulfanilamide,
#'             urethane, temperature and pressure. Journal of cellular and comparative
#'             physiology 20, 247-268.
#' 
#'             Kattge J., Knorr W., 2007: Temperature acclimation in a biochemical
#'             model of photosynthesis: a reanalysis of data from 36 species.
#'             Plant, Cell and Environment 30, 1176-1190.
Arrhenius.temp.response <- function(param,Temp,Ha=NULL,Hd=NULL,dS=NULL,
                                    constants=bigleaf.constants()){
  
  Temp <- Temp + constants$Kelvin
  Tref <- 25.0 + constants$Kelvin
  
  Ha <- Ha * 1000
  Hd <- Hd * 1000
  
  if (is.null(Ha)){
    
    stop("Activation energy (Ha) has to be provided!")
    
  }
  
  if (is.null(Hd) & is.null(dS)){
    
    param25 <- param / exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
  
  } else if (!is.null(Hd) & !is.null(dS)){
    
    param25 <- param / 
      ( exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp)) *
          (1 + exp((Tref*dS - Hd) / (Tref * constants$Rgas))) /
          (1 + exp((Temp*dS - Hd) / (Temp * constants$Rgas)))
      )
    
  } else if ((!is.null(Hd) & is.null(dS)) | (is.null(Hd) & !is.null(dS)) ){
    
    warning("Both Hd and dS have to be provided for a temperature response
             that considers a temperature optimum and a deactivation term!
             Continue considering activation energy (Ha) only...")
    
    param25 <- param / exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
    
  }

  return(param25)
  
}



#' Stomatal slope parameter "g1"
#' 
#' @description Estimation of the intrinsic WUE metric "g1" (stomatal slope) 
#'              from nonlinear regression.
#' 
#' @param data       Data.frame or matrix containing all required columns
#' @param Tair       Air temperature (deg C)
#' @param pressure   Atmospheric pressure (kPa)
#' @param GPP        Gross primary productivity (umol CO2 m-2 s-1)
#' @param Gs         Surface conductance to water vapor (mol m-2 s-1)
#' @param VPD        Vapor pressure deficit (kPa)
#' @param Ca         Atmospheric CO2 concentration (umol mol-1)
#' @param model      Stomatal model used. One of c("USO","Ball&Berry","Leuning")
#' @param robust.nls Use robust nonlinear regression (\code{\link[robustbase]{nlrob}})? Default is FALSE.
#' @param nmin       Minimum number of data required to perform the fit; defaults to 40.
#' @param fitg0      Should g0 and g1 be fitted simultaneously? 
#' @param g0         Minimum stomatal conductance (mol m-2 s-1); ignored if \code{fitg0} is TRUE.
#' @param fitD0      Should D0 be fitted along with g1 (and g0 if fitg0 = TRUE)?; only used if \code{model} is "Leuning"
#' @param D0         Stomatal sensitivity parameter to VPD; only used if \code{model} is "Leuning" and fitD0 is FALSE
#' @param Gamma      Canopy CO2 compensation point (umol mol-1); only used if \code{model} is "Leuning". 
#'                   Can be a constant or a variable. Defaults to 50 umol mol-1.
#' @param constants  Kelvin - conversion degree Celsius to Kelvin \cr
#'                   Rgas - universal gas constant (J mol-1 K-1)
#' 
#' @details All stomatal models were developed at leaf-level, but its parameters 
#'          can also be estimated at ecosystem level (but be aware of caveats).
#'          
#'          The unified stomatal optimization (USO) model is given as (Medlyn et al. 2011):
#'      
#'          \deqn{gs = g0 + 1.6*(1.0 + g1/sqrt(VPD)) * GPP/Ca}
#'          
#'          The semi-empirical model by Ball et al. 1987 is defined as:
#'          
#'          \deqn{gs = g0 + g1* ((An * rH) / Ca)}
#'          
#'          Leuning 1995 suggested a revised version of the Ball&Berry model:
#'          
#'          \deqn{gs = g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0))}
#'          
#'          The parameters in the models are estimated using nonlinear regression (\code{\link[stats]{nls}}) if
#'          \code{robust.nls == FALSE} and weighted nonlinear regression if \code{robust.nls == TRUE}.
#'          The weights are calculated from \code{\link[robustbase]{nlrob}}, and \code{\link[stats]{nls}}
#'          is used for the actual fitting.
#'          Alternatively to measured VPD and Ca (i.e. conditions at instrument height), conditions at 
#'          the big-leaf surface can be provided. They can be calculated using \code{\link{surface.conditions}}.
#'          
#' 
#' @return an nls model object, containing information on the fitted parameters, their uncertainty range,
#'         model fit, etc.
#' 
#' @references Medlyn B.E., et al., 2011: Reconciling the optimal and empirical approaches to
#'             modelling stomatal conductance. Global Change Biology 17, 2134-2144.
#'             
#'             Ball T.J., Woodrow I.E., Berry J.A. 1987: A model predicting stomatal conductance
#'             and its contribution to the control of photosynthesis under different environmental conditions.
#'             In: Progress in Photosynthesis Research, edited by J.Biggins, pp. 221-224, Martinus Nijhoff Publishers,
#'             Dordrecht, Netherlands.
#'             
#'             Leuning R., 1995: A critical appraisal of a combined stomatal-photosynthesis
#'             model for C3 plants. Plant, Cell and Environment 18, 339-355.
#' 
#' @examples 
#' ## filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_Jun_2014_2 <- filter.data(DE_Tha_Jun_2014,quality.control=FALSE,
#'                                  vars.qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter.growseas=FALSE,filter.precip=TRUE,
#'                                  filter.vars=c("Tair","PPFD","ustar","LE"),
#'                                  filter.vals.min=c(5,200,0.2,0),
#'                                  filter.vals.max=c(NA,NA,NA,NA),NA.as.invalid=TRUE,
#'                                  quality.ext="_qc",good.quality=c(0,1),
#'                                  missing.qc.as.bad=TRUE,GPP="GPP_nt",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
#'                                  tprecip=0.1,precip.hours=24,records.per.hour=2)
#'                                  
#' DE_Tha_Jun_2014_2[DE_Tha_Jun_2014_2[,"valid"] < 1,] <- NA
#' 
#' # calculate Gs from the the inverted PM equation
#' Ga <- aerodynamic.conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default).
#' Gs_PM <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,PM=TRUE)[,"Gs_mol"]
#'                              
#' ### Estimate the stomatal slope parameter g1 using the USO model
#' mod_USO <- stomatal.slope(DE_Tha_Jun_2014_2,model="USO",GPP="GPP_nt",Gs=Gs_PM,
#'                           robust.nls=FALSE,nmin=40,fitg0=FALSE)
#'                           
#' ### Use robust regression to minimize influence of outliers in Gs                           
#' mod_USO <- stomatal.slope(DE_Tha_Jun_2014_2,model="USO",GPP="GPP_nt",Gs=Gs_PM,
#'                           robust.nls=TRUE,nmin=40,fitg0=FALSE)
#' 
#' ### Estimate the same parameter from the Ball&Berry model and prescribe g0
#' mod_BB <- stomatal.slope(DE_Tha_Jun_2014_2,model="Ball&Berry",GPP="GPP_nt",
#'                          robust.nls=FALSE,Gs=Gs_PM,g0=0.01,nmin=40,fitg0=FALSE)
#' 
#' ## same for the Leuning model, but this time estimate both g1 and g0 (but fix D0)
#' mod_Leu <- stomatal.slope(DE_Tha_Jun_2014_2,model="Leuning",GPP="GPP_nt",Gs=Gs_PM,
#'                           robust.nls=FALSE,nmin=40,fitg0=FALSE,D0=1.5,fitD0=FALSE)
#' 
#' @importFrom stats nls na.exclude
#' @importFrom robustbase nlrob 
#' 
#' @export 
stomatal.slope <- function(data,Tair="Tair",pressure="pressure",GPP="GPP",Gs="Gs",
                           VPD="VPD",Ca="Ca",model=c("USO","Ball&Berry","Leuning"),
                           robust.nls=FALSE,nmin=40,fitg0=FALSE,g0=0,fitD0=FALSE,
                           D0=1.5,Gamma=50,constants=bigleaf.constants()){
  
  model <- match.arg(model)
  
  check.input(data,list(Tair,pressure,GPP,Gs,VPD,Ca))
  
  df <- data.frame(Tair,pressure,GPP,Gs,VPD,Ca)
  
  if (model == "Leuning"){
    if (length(Gamma) == 1){
      Gamma <- rep(Gamma,nrow(df))
    }
    df$Gamma <- Gamma
  }
  
  nr_data <- sum(!is.na(GPP) & !is.na(Gs) & !is.na(VPD) & !is.na(Ca))
  
  if (nr_data < nmin){
    stop("number of data is less than 'nmin'. g1 is not fitted to the data.")
  } else {
    
    if (model == "USO"){
      
      if (fitg0){
        if (robust.nls){
          mod_weights <- nlrob(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g0=0,g1=3),
                               na.action=na.exclude)$w
          mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3),weights=mod_weights)
        } else {
          mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3))
        }
      } else {
        if (robust.nls){
          df$g0   <- rep(g0,nrow(df)) # g0 as constant does not work in the nlrob function...
          mod_weights <- nlrob(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g1=3),
                               na.action=na.exclude)$w
          mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3),weights=mod_weights)
        } else {
          mod <- nls(Gs ~ g0 + 1.6*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3))
        }
      }
      
    } else if (model == "Leuning"){
      
      if (fitg0){
        if (fitD0){
          if (robust.nls){
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9,D0=1.5),na.action=na.exclude)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5),
                       weights=mod_weights)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5))
          }
        } else {
          if (robust.nls){
            df$D0  <- rep(D0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9),na.action=na.exclude)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9),
                       weights=mod_weights)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9))
          }
        }
      } else {
        if (fitD0){
          if (robust.nls){
            df$g0   <- rep(g0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9,D0=1.5),na.action=na.exclude)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5),
                       weights=mod_weights)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5))
          }
        } else {
          if (robust.nls){
            df$g0  <- rep(g0,nrow(df))
            df$D0  <- rep(D0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9),na.action=na.exclude)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9),
                       weights=mod_weights)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9))
          }
        }
      }
      
    } else if (model == "Ball&Berry"){
      
      rH <- VPD.to.rH(VPD,Tair)
      df$rH <- rH
      
      if (fitg0){
        if (robust.nls){
          mod_weights <- nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),data=df,
                               na.action=na.exclude)$w
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),weights=mod_weights)
        } else {
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9))
        }
      } else {
        if (robust.nls){
          g0   <- rep(g0,nrow(df))
          mod_weights <- nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),data=df,
                               na.action=na.exclude)$w
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),weights=mod_weights)
        } else {
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9))
        }
      }
      
    }
    
  }
  
  return(mod)
  
}






#' Ecosytem light response
#' 
#' @description calculates GPPmax at a reference (mostly saturating) PPFD and 
#'              ecosystem quantum yield using a rectangular light response curve.
#' 
#' @param data      Data.frame or matrix containing all required columns
#' @param NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param Reco      Ecosystem respiration (umol CO2 m-2 s-1)
#' @param PPFD      Photosynthetic photon flux density (umol m-2 s-1)
#' @param PPFD_ref  Reference PPFD (umol m-2 s-1) for which GPPmax is estimated.
#'                  Default is 2000 umol m-2 s-1.
#' 
#' @details A rectangular light response curve is fitted to NEE data. The curve
#'          takes the form as described in Falge et al. 2001:
#'          
#'          \deqn{NEE = \alpha PPFD / (1 - (PPFD / PPFD_ref) + \alpha 
#'                       PPFD / GPPmax)- Reco}
#'                       
#'          where \code{\alpha} is the ecosystem quantum yield (umol CO2 m-2 s-1) (umol quanta m-2 s-1)-1, 
#'          and GPPmax is the GPP at the reference PPFD (usually at saturating light). \code{\alpha} 
#'          represents the slope of the light response curve, and is a measure for the light use
#'          efficiency of the canopy. 
#'          
#'          The advantage of this equation over the standard rectangular light response
#'          curve is that GPPmax at PPFD_ref is more readily interpretable
#'          as it constitutes a value observed in the ecosystem, in contrast to 
#'          GPPmax (mostly named 'beta') in the standard model that occurs at infinite light.
#'          PPFD_ref defaults to 2000 umol m-2 s-1, but other values can be used. For 
#'          further details refer to Falge et al. 2001.
#' 
#' @note   Note the sign convention. Negative NEE indicates that carbon is taken up
#'         by the ecosystem. Reco has to be 0 or larger.
#' 
#' @return A nls model object containing estimates (+/- SE) for alpha and GPPmax.
#' 
#' @references Falge E., et al. 2001: Gap filling strategies for defensible annual
#'             sums of net ecosystem exchange. Agricultural and Forest Meteorology 107,
#'             43-69.
#' 
#' @importFrom stats nls
#' 
#' @export
light.response <- function(data,NEE,Reco,PPFD,PPFD_ref=2000){

  check.input(data,list(NEE,Reco,PPFD))
  
  mod <- nls(-NEE ~ alpha * PPFD / (1 - (PPFD / PPFD_ref) + (alpha * PPFD / GPP_max)) + Reco,
             start=list(alpha=0.02,GPP_max=30))
  
  return(mod)
}  


  

#' Light use efficiency (LUE)
#' 
#' @description amount of carbon fixed (GPP) per incoming light.
#' 
#' @param GPP     Gross ecosystem productivity (umol CO2 m-2 s-1)
#' @param PPFD    Photosynthetic photon flux density (umol quanta m-2 s-1)
#' 
#' @details Light use efficiency is calculated as
#'          
#'          \deqn{LUE = sum(GPP)/sum(PPFD)}
#'          
#'          where both GPP and PPFD are in umol m-2 s-1. A more meaningful 
#'          (as directly comparable across ecosystems) is to take absorbed PPFD
#'          rather than incoming PPFD as used here.
#' 
#' @return \item{LUE -}{Light use efficiency (-)}
#' 
#' @examples
#' light.use.efficiency(GPP=20,PPFD=1500)
#' 
#' @export 
light.use.efficiency <- function(GPP,PPFD){
  
  comp <- complete.cases(GPP,PPFD)
  
  LUE <- sum(GPP[comp],na.rm=TRUE)/sum(PPFD[comp],na.rm=TRUE)
  
  return(c("LUE"=LUE))
}
  



  
  
  