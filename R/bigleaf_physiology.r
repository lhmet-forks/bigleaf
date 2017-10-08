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
#' @param RecoLeaf   Ecosytem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0          
#' @param Gs         Surface conductance to water vapor (mol m-2 s-1)
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
intercellular.CO2 <- function(data,Ca,GPP,RecoLeaf=NULL,Gs,calc.Csurf=FALSE,
                              Ga_CO2=NULL,NEE=NULL,Tair=NULL,pressure=NULL,
                              constants=bigleaf.constants()){
  
  Ca  <- check.columns(data,Ca)
  GPP <- check.columns(data,GPP)
  Gs  <- check.columns(data,Gs)
  
  if (calc.Csurf){
    Tair     <- check.columns(data,Tair)
    pressure <- check.columns(data,pressure)
    Ga_CO2   <- check.columns(data,Ga_CO2)
    NEE      <- check.columns(data,NEE)
    
    Ca <- Ca.surface(Ca,NEE,Ga_CO2,Tair,pressure)
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
#' @param PPFD_c    PPFD threshold, above which the canopy is considered to 
#'                  be light saturated (and Rubisco limited).
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
#'         \item{Vcmax}{maximum carboxylation rate at Temp (umol m-2 s-1)}
#'         \item{Vcmax25}{maximum carboxylation rate at 25degC (umol m-2 s-1)}
#'         \item{Ko}{Michaelis-Menten constant for O2 at Temp (mmol mol-1)}
#'         \item{Kc}{Michaelis-Menten constant for CO2 at Temp (umol mol-1)}   
#'        
#' @references Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385. 
#'             
#'             Bernacchi C.J, Singsaas E.L., Pimentel C., Portis JR A.R., Long S.P., 2001:
#'             Improved temperature response functions for models of Rubisco-limited
#'             photosynthesis. Plant, Cell and Environment 24, 253-259. 
#'             
#'             Kattge J., Knorr W., 2007: Temperature acclimation in a biochemical
#'             model of photosynthesis: a reanalysis of data from 36 species.
#'             Plant, Cell and Environment 30, 1176-1190.
#'             
#'
#' @examples 
#' carbox.rate(Temp=20,GPP=15,Ci=300,PPFD=1000,PPFD_c=500)  
#'                                       
#' @export                  
carbox.rate <- function(Temp,GPP,Ci,PPFD,PPFD_c,Oi=0.21,Kc25=404.9,Ko25=278.4,
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
  GPP[PPFD < PPFD_c] <- NA
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
stomatal.slope <- function(data,Tair="Tair",pressure="pressure",GPP="GPP_nt",Gs="Gs",
                           VPD="VPD",Ca="Ca",model=c("USO","Ball&Berry","Leuning"),
                           robust.nls=FALSE,nmin=40,fitg0=FALSE,g0=0,fitD0=FALSE,
                           D0=1.5,Gamma=50,constants=bigleaf.constants()){
  
  model <- match.arg(model)
  
  Tair     <- check.columns(data,Tair)
  pressure <- check.columns(data,pressure)
  GPP      <- check.columns(data,GPP)
  Gs       <- check.columns(data,Gs)
  VPD      <- check.columns(data,VPD)
  Ca       <- check.columns(data,Ca)
  check.length(Tair,pressure,GPP,Gs,VPD,Ca)
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