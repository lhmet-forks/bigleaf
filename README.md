# bigleaf #

bigleaf is an R package for the calculation of of physical (e.g. aerodynamic conductance, surface temperature) and physiological
(e.g. canopy conductance, water-use efficiency) ecosystem properties from eddy covariance data and accompanying meteorological measurements. 
All calculations are based on a 'big-leaf' representation of the vegetation, in which vertical meteorological variations within the canopy
are ignored and all fluxes are assumed to originate from a single horizontal plane within the canopy.

Installation: the package is on CRAN and can be installed using:
install.packages("bigleaf")

Example:
Calculate surface conductance for a spruce forest
surface.conductance(DE_Tha_June_2014,Tair="Tair",pressure="pressure",Rn="Rn",VPD="VPD",LE="LE",Ga="Ga")
surface.conductance(DE_Tha_June_2014,Tair="Tair",pressure="pressure",Rn="Rn",VPD="VPD",LE="LE",Ga=0.1)

Here, DE_Tha_June_2014 denotes the input data.frame/matrix. Note that
input variables can be provided as column names of the input data.frame/matrix (as argument Ga in line 1 above),
or alternatively, as vectors with the same length as the input data.frame/matrix or of length 1 (as argument Ga in line 2 above).

The package provides the following functionalities:

# data filtering
- data quality filter
- filter based on meteorological variables (radiation, precipitation, ustar, temperature, etc.)
- growing season filter (based on daily GPP)

# meteorological variables:
- air density
- virtual temperature
- pressure from altitude
- psychrometric constant
- latent heat of vaporization
- saturation vapor pressure
- slope of saturation vapor pressure curve

# aerodynamic properties:
- aerodynamic conductance (different versions) for water, heat, momentum, and CO2
- Canopy boundary-layer conductance (Rb and kB-1 parameter; different models)
- Monin-Obhukov length
- stability parameter zeta
- stability correction functions (different versions)
- roughness length for momentum (z0m) estimation
- decoupling coefficient 'omega'
- Reynolds number
- wind speed at given height from wind profile equation

# Evapotranspiration (ET) and water-use efficiency (WUE)
- potential ET (Priestley-Taylor equation)
- reference ET (Penman-Monteith equation)
- imposed and equilibrum ET
- WUE, inherent WUE, underlying WUE

# physiological variables:
- canopy conductance (inverted Penman-Monteith equation)
- bulk intercellulary CO2 concentrtation (Ci)
- bulk photosynthetic capacity (Vcmax25 and Jmax25)
- stomatal slope g1 (USO, Ball&Berry, Leuning models)
- stomatal sensitivity to VPD
- ecosystem light response, light-use efficiency

# surface conditions
- VPD, Ca, Temperature, vapor pressure, specific humidity at the big-leaf surface

# energy balance
- biochemical energy
- energy balance ratio (EBR)
- "missing" energy

# unit conversions
- conductance conversion from ms-1 to mol m-2 s-1
- conversions between humidity measures vapor pressure, specific humidity, relative humidity, and VPD
- conversion between laten heat flux and evapotranspiration