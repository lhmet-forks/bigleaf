################################
### Aerodynamic conductance #### 
################################

#' Aerodynamic Conductance
#' 
#' @description Bulk aerodynamic conductance, including options for the boundary layer conductance
#'              formulation and stability correction functions.
#' 
#' @param data              Data.frame or matrix containing all required variables
#' @param Tair              Air temperature (deg C)
#' @param pressure          Atmospheric pressure (kPa)
#' @param wind              Wind speed (m s-1)
#' @param ustar             Friction velocity (m s-1)
#' @param H                 Sensible heat flux (W m-2)
#' @param zr                Instrument (reference) height (m)
#' @param zh                Canopy height (m)
#' @param d                 Zero-plane displacement height (m)
#' @param z0m               Roughness length for momentum (m)
#' @param Dl                Characteristic leaf dimension (m); only used if \code{Rb_model} is \code{"Choudhury_1988"} or \code{"Su_2001"}.
#' @param N                 Number of leaf sides participating in heat exchange (1 or 2); only used if \code{Rb_model = "Su_2001"}.
#'                          Defaults to 2.
#' @param fc                Fractional vegetation cover (-); only used if \code{Rb_model = "Su_2001"}. See Details.
#' @param LAI               One-sided leaf area index (m2 m-2); only used if \code{Rb_model} is \code{"Choudhury_1988"} or \code{"Su_2001"}.
#' @param Cd                Foliage drag coefficient (-); only used if \code{Rb_model = "Su_2001"}. 
#' @param hs                Roughness length of bare soil (m); only used if \code{Rb_model = "Su_2001"}.
#' @param wind_profile      Should Ga for resistance be calculated based on the logarithmic wind profile equation? 
#'                          Defaults to \code{FALSE}.
#' @param stab_correction   Should stability correction be applied? Defaults to \code{TRUE}. Ignored if \code{wind_profile = FALSE}.                         
#' @param stab_formulation  Stability correction function. Either \code{"Dyer_1970"} (default) or
#'                          \code{"Businger_1971"}. Ignored if \code{wind_profile = FALSE} or if \code{stab_correction = FALSE}.
#' @param Rb_model          Boundary layer resistance formulation. One of \code{"Thom_1972","Choudhury_1988","Su_2001","constant_kB-1"}.
#' @param kB                kB-1 value; only used if \code{Rb_model = "constant_kB-1"}
#' @param Sc                Optional: Schmidt number of additional quantities to be calculated
#' @param Sc_name           Optional: Name of the additonal quantities, has to be of same length than 
#'                          \code{Sc_name}
#' @param constants         k - von Karman constant \cr
#'                          cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                          Kelvin - conversion degree Celsius to Kelvin \cr
#'                          g - gravitational acceleration (m s-2) \cr
#'                          pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                          Tair0 - reference air temperature (K) \cr
#'                          Sc_CO2 - Schmidt number for CO2 \cr 
#'                          Pr - Prandtl number (if \code{Sc} is provided)
#'
#' @details
#' 
#' Aerodynamic conductance for heat (Ga_h) is calculated as:
#' 
#'     \deqn{Ga_h = 1 / (Ra_m + Rb)}
#'  
#'  where Ra_m is the aerodynamic resistance for momentum and Rb the (quasi-laminar) canopy boundary
#'  layer resistance ('excess resistance').
#'  
#'  The aerodynamic resistance for momentum Ra_m is given by:
#'  
#'     \deqn{Ra_m = u/ustar^2}
#'  
#'  Note that this formulation accounts for changes in atmospheric stability, and does not require an
#'  additional stability correction function. 
#'  
#'  An alternative method to calculate Ra_m is provided
#'  (calculated if \code{wind_profile = TRUE}):
#'  
#'     \deqn{Ra_m = (ln((zr - d)/z0m) - psi_h) / (k ustar)}
#'  
#'  If the roughness parameters z0m and d are unknown, they can be estimated using
#'  \code{\link{roughness.parameters}}. The argument \code{stab_formulation} determines the stability correction function used 
#'  to account for the effect of atmospheric stability on Ra_m (Ra_m is lower for unstable
#'  and higher for stable stratification). Stratification is based on a stability parameter zeta (z-d/L),
#'  where z = reference height, d the zero-plane displacement height, and L the Monin-Obukhov length, 
#'  calculated with \code{\link{Monin.Obukhov.length}}
#'  The stability correction function is chosen by the argument \code{stab_formulation}. Options are 
#'  \code{"Dyer_1970"} and \code{"Businger_1971"}.
#'  
#'  The model used to determine the canopy boundary layer resistance (Rb) is specified by 
#'  the argument \code{Rb_model}. The following options are implemented:
#'  \code{"Thom_1972"} is an empirical formulation based on the friction velocity (ustar) (Thom 1972):
#'  
#'    \deqn{Rb = 6.2ustar^-0.667}
#'    
#'  The model by Choudhury 1988 (\code{"Choudhury_1988"}), calculates Rb
#'  based on leaf width, LAI and ustar (Note that in the original formulation leaf width
#'  instead of the characteristic leaf dimension (Dl) is taken):
#'   
#'     \deqn{Gb = LAI((2a/\alpha)*sqrt(u(h)/w)*(1-exp(-\alpha/2)))}
#'     
#'  The option \code{"Su_2001"} calculates Rb based on the physically-based Rb model by Su et al. 2001,
#'  a simplification of the model developed by Massman 1999:
#'  
#'     \deqn{kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2}
#'  
#'  The models calculate the parameter kB-1, which is related to Rb:
#'  
#'     \deqn{kB-1 = Rb * (k * ustar)}
#'  
#'  Both kB-1 and Rb, as well as Gb are given in the output.
#'  
#'  Rb (and Gb) for water vapor and heat are assumed to be equal. Rb for other quantities
#'  x is calculated as (Hicks et al. 1987):
#'  
#'     \deqn{Rb_x = Rb * (Sc_x / Pr)^0.67}
#'  
#'  where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
#'  
#'  data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Gb,Rb,kB,zeta,psi_h,Ra_x,Rb_x))
#' @return a dataframe with the following columns:
#'         \item{Ga_m}{Aerodynamic conductance for momentum (m s-1)}
#'         \item{Ra_m}{Aerodynamic resistance for momentum (s m-1)}
#'         \item{Ga_h}{Aerodynamic conductance for heat and water vapor (m s-1)}
#'         \item{Ra_h}{Aerodynamic resistance for heat and water vapor (s m-1)}
#'         \item{Gb}{Canopy boundary layer conductance for heat and water vapor (m s-1)}
#'         \item{Rb}{Canopy boundary layer resistance for heat and water vapor (s m-1)}
#'         \item{kB}{kB-1 parameter}
#'         \item{zeta}{Stability parameter 'zeta'}
#'         \item{psi_h}{Integrated stability correction function}
#'         \item{Ga_CO2}{Aerodynamic conductance for CO2 (m s-1)}
#'         \item{Ra_CO2}{Aerodynamic resistance for CO2 (s m-1)}
#'         \item{Ra_Sc_name}{Aerodynamic resistance of \code{Sc_name} (s m-1). Only added if \code{Sc_name} and 
#'                           \code{Sc_name} are provided}
#'         \item{Rb_CO2}{Canopy boundary layer resistance for CO2 (s m-1)}
#'         \item{Rb_Sc_name}{Boundary layer resistance of \code{Sc_name} (s m-1). Only added if \code{Sc_name} and 
#'                           \code{Sc_name} are provided}
#'        
#' @note The roughness length for water and heat (z0h) is not returned by this function, but 
#'       it can be calculated from the following relationship (e.g. Verma 1989):
#'       
#'       \deqn{kB-1 = ln(z0m/z0h)} 
#'       
#'       it follows:
#'       
#'       \deqn{z0h = z0m / exp(kB-1)}
#'       
#'       \code{kB-1} is an output of this function.
#'       
#'       Note that input variables such as LAI, Dl, or zh can be either constants, or
#'       time varying (i.e. vectors of the same length as \code{data}).
#'         
#' @references Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
#'             In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
#'             
#'             Verhoef, A., De Bruin, H., Van Den Hurk, B., 1997: Some practical notes on the parameter kB-1
#'             for sparse vegetation. Journal of Applied Meteorology, 36, 560-572.
#'             
#'             Hicks, B.B., Baldocchi, D.D., Meyers, T.P., Hosker, J.R., Matt, D.R., 1987:
#'             A preliminary multiple resistance routine for deriving dry deposition velocities
#'             from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
#'             
#'             Monteith, J.L., Unsworth, M.H., 2008: Principles of environmental physics.
#'             Third Edition. Elsevier Academic Press, Burlington, USA. 
#' 
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Choudhury}}, \code{\link{Gb.Su}} for calculations of Rb / Gb only
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
                                    zr,zh,d,z0m,Dl,N=2,fc=NULL,LAI,Cd=0.2,hs=0.01,wind_profile=FALSE,
                                    stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                    Rb_model=c("Thom_1972","Choudhury_1988","Su_2001","constant_kB-1"),
                                    kB=NULL,Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){
  
  Rb_model         <- match.arg(Rb_model)
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,list(Tair,pressure,wind,ustar,H))

  ## calculate canopy boundary layer conductance (Gb)
  if (Rb_model %in% c("Thom_1972","Choudhury_1988","Su_2001")){
    
    if (Rb_model == "Thom_1972"){
      
      Gb_mod <- Gb.Thom(ustar=ustar,Sc=Sc,Sc_name=Sc_name,constants=constants)
      
    } else if (Rb_model == "Choudhury_1988"){
      
      Gb_mod <- Gb.Choudhury(data,Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,
                             H=H,leafwidth=Dl,LAI=LAI,zh=zh,zr=zr,d=d,
                             stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                             constants=constants)
      
    } else if (Rb_model == "Su_2001"){
      
      Gb_mod <- Gb.Su(data=data,Tair=Tair,pressure=pressure,ustar=ustar,wind=wind,
                      H=H,zh=zh,zr=zr,d=d,Dl=Dl,N=N,fc=fc,LAI=LAI,Cd=Cd,hs=hs,
                      stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                      constants=constants)  
      
    }
    
    kB   <- Gb_mod[,"kB"]
    Rb   <- Gb_mod[,"Rb"]
    Rb_x <- data.frame(Gb_mod[,grep(colnames(Gb_mod),pattern="Rb_")])
    colnames(Rb_x) <- grep(colnames(Gb_mod),pattern="Rb_",value=TRUE)
    Gb   <- Gb_mod[,"Gb"]
    
  } else if (Rb_model == "constant_kB-1"){
    
    if(is.null(kB)){
      stop("value of kB-1 has to be specified if Rb_model is set to 'constant_kB-1'!")
    } else {
      Rb     <- kB / (constants$k * ustar)
      Gb     <- 1 / Rb
      
      if (!is.null(Sc) | !is.null(Sc_name)){
        if (length(Sc) != length(Sc_name)){
          stop("arguments 'Sc' and 'Sc_name' must have the same length")
        }
        if (!is.numeric(Sc)){
          stop("argument 'Sc' must be numeric")
        }
      }
      
      Sc   <- c(constants$Sc_CO2,Sc)
      Rb_x <- data.frame(lapply(Sc,function(x) Rb * (x/constants$Pr)^0.67))
      colnames(Rb_x) <- paste0("Rb_",c("CO2",Sc_name))
      
    }
    
  } 
  
  ## calculate aerodynamic conductance for momentum (Ga_m)
  if (wind_profile){
    
    if (stab_correction){
      
      zeta  <-  stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,
                                    H=H,zr=zr,d=d,constants=constants)
      
      if (stab_formulation %in% c("Dyer_1970","Businger_1971")){
        
        psi_h <- stability.correction(zeta,formulation=stab_formulation)[,"psi_h"]
        
      } else {
        stop("'stab_formulation' has to be one of 'Dyer_1970' or 'Businger_1971'.
             Choose 'stab_correction = FALSE' if no stability correction should be applied.")
      }
      
      Ra_m  <- pmax((log((zr - d)/z0m) - psi_h),0) / (constants$k*ustar)
      
      } else {
        
        Ra_m  <- pmax((log((zr - d)/z0m)),0) / (constants$k*ustar)
        zeta = psi_h <- rep(NA_integer_,length=length(Ra_m))
        
      }
    
  } else {
    
    if ((!missing(zr) | !missing(d) | !missing(z0m)) & Rb_model %in% c("constant_kB-1","Thom_1972")){
      warning("Provided roughness length parameters (zr,d,z0m) are not used if 'wind_profile = FALSE' (the default). Ra_m is calculated as wind / ustar^2")
    }
    
    Ra_m <- wind / ustar^2
    zeta = psi_h <- rep(NA_integer_,length=length(Ra_m))
    
  }
  
  Ga_m   <- 1/Ra_m
  Ra_h   <- Ra_m + Rb
  Ga_h   <- 1/Ra_h
  Ra_x   <- Ra_m + Rb_x
  Ga_CO2 <- 1/Ra_x[,1]
  colnames(Ra_x) <- paste0("Ra_",c("CO2",Sc_name))

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Gb,Rb,kB,zeta,psi_h,Ga_CO2,Ra_x,Rb_x))
  
}