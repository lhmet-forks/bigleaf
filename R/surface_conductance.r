##############################
#### Surface conductance  ####
##############################

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
#' @details Surface conductance (Gs) is calculated from the inverted Penman-Monteith equation 
#'          (if \code{PM} is TRUE, the default):
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
#' # calculate Gs based on a simple gradient approach
#' Gs_gradient <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                                    VPD="VPD",PM=FALSE)
#' summary(Gs_gradient)
#' 
#' # calculate Gs from the the inverted PM equation (now Rn, and Ga are needed),
#' # using a simple estimate of Ga based on Thom 1972
#' Ga <- aerodynamic.conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default).
#' # Note that Ga is not added to the data.frame 'DE_Tha_Jun_2014'
#' Gs_PM <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,PM=TRUE)
#' summary(Gs_PM)
#' 
#'                               
#' # now add Ga to the data.frame 'DE_Tha_Jun_2014' and repeat
#' DE_Tha_Jun_2014_2$Ga <- Ga
#' Gs_PM2 <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                               Rn="Rn",G="G",S=NULL,VPD="VPD",Ga="Ga",PM=TRUE)
#' # note the difference to the previous version (Ga="Ga")
#' summary(Gs_PM2)
#'                              
#' @export
surface.conductance <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga",missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                PM=TRUE,constants=bigleaf.constants()){ 
  
  check.input(data,list(Tair,pressure,VPD,LE))
  
  if (!PM){
    
    Gs_simple_mol <- (LE.to.ET(LE,Tair)/constants$Mw) * pressure / VPD
    Gs_simple_ms  <- mol.to.ms(Gs_simple_mol,Tair,pressure)
    
    return(data.frame(Gs_ms=Gs_simple_ms,Gs_mol=Gs_simple_mol))
    
  } else {
    
    check.input(data,list(Rn,Ga))

    if(!is.null(G)){
      check.input(data,list(G))
      if (!missing.G.as.NA){G[is.na(G)] <- 0}
    } else {
      warning("Ground heat flux G is not provided and set to 0.")
      G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
    }
    
    if(!is.null(S)){
      check.input(data,list(S))
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