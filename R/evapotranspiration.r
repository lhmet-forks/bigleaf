###########################
#### Evapotranspiration ###
###########################

#' Potential evapotranspiration
#' 
#' @description Potential evapotranspiration ET_pot according to Priestley & Taylor 1972.
#' 
#' @param data      Data.frame or matrix containing all required variables; optional
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
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
                   missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,Rn))
  
  if(!is.null(G)){
    check.input(data,list(G))
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.",call.=FALSE)
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    check.input(data,list(S))
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    warning("Energy storage fluxes S are not provided and set to 0.",call.=FALSE)
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
#' @param pressure  Atmospheric pressure (kPa)
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
#' # Calculate ET_ref for a surface with known Gs (0.5mol m-2 s-1) and Ga (0.1 ms-1)
#' 
#' # Gs is required in ms-1
#' Gs_ms <- mol.to.ms(0.5,Tair=20,pressure=100)
#' ET_ref <- ET.ref(Gs=Gs_ms,Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400)
#' 
#' # now cross-check with the inverted version
#' surface.conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=ET_ref[,"LE_ref"])
#' 
#' @export                 
ET.ref <- function(data,Gs=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga",
                   G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                   constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,VPD,Rn,Ga))
  
  if(!is.null(G)){
    check.input(data,list(G))
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    check.input(data,list(S))
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
#' @param pressure  Atmospheric pressure (kPa)
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
#'       Aerodynamic conductance (Ga) can be calculated using \code{\link{aerodynamic.conductance}}.
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
#' @examples 
#' df <- data.frame(Tair=20,pressure=100,VPD=seq(0.5,4,0.5),
#'                  Gs=seq(0.01,0.002,length.out=8),Rn=seq(50,400,50))            
#' ET.components(df)            
#'             
#' @export
ET.components <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs",
                          Rn="Rn",G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                          constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,VPD,Rn,Gs))
  
  if(!is.null(G)){
    check.input(data,list(G))
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    warning("ground heat flux G is not provided and set to 0.")
    G <- rep(0,ifelse(!missing(data),nrow(data),length(Tair)))
  }
  
  if(!is.null(S)){
    check.input(data,list(S))
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