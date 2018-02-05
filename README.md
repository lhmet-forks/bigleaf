# The bigleaf R package #

![](https://www.r-pkg.org/badges/version/bigleaf) ![](http://cranlogs.r-pkg.org/badges/grand-total/bigleaf)

**bigleaf** is an R package for the calculation of physical (e.g. aerodynamic conductance, surface temperature) and physiological
(e.g. canopy conductance, water-use efficiency) ecosystem properties from eddy covariance data and accompanying meteorological measurements. 
All calculations are based on a 'big-leaf' representation of the vegetation and return representative bulk ecosystem/canopy variables.

# Installation
The bigleaf R package is on CRAN and can be installed using:
```
install.packages("bigleaf")
```

The development version can be directly installed from this bitbucket repository: 
```
library(devtools)
install_bitbucket("juergenknauer/bigleaf")
```


# Usage
Most functions work by providing a data.frame or matrix which contains all required variables:
For example, surface conductance for the spruce forest in Tharandt, Germany (DE-Tha) can be calculated with 
the following commands:
```
DE_Tha_June_2014$Ga <- aerodynamic.conductance(DE_Tha_June_2014,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar")[,"Ga_h"]
surface.conductance(DE_Tha_June_2014,Tair="Tair",pressure="pressure",Rn="Rn",VPD="VPD",LE="LE",Ga="Ga")
surface.conductance(DE_Tha_June_2014,Tair="Tair",pressure="pressure",Rn="Rn",VPD="VPD",LE="LE",Ga=0.1)
```
Here, DE_Tha_June_2014 denotes the input data.frame. Note that input variables can be provided as column names of the 
input data.frame (as argument Ga in line 2 above), or alternatively, as vectors with the same length as the input data.frame
or of length 1 (as argument Ga in line 3 above). If variables are provided in the default column names (as above), the command can 
be shortened to:
```
surface.conductance(DE_Tha_June_2014)
```
Important: please ensure that all input variables are in the correct units as described on the help pages.

[Please report bugs or issues here](https://bitbucket.org/juergenknauer/bigleaf/issues?status=new&status=open)

# Package content 
The package provides the following functionalities:

## Data filtering
- data quality filter
- filter based on meteorological variables (radiation, precipitation, ustar, temperature, etc.)
- growing season filter (based on daily GPP)

## Meteorological variables
- air density
- virtual temperature
- pressure from altitude
- psychrometric constant
- latent heat of vaporization
- saturation vapor pressure
- slope of saturation vapor pressure curve

## Aerodynamic properties
- aerodynamic conductance (different versions) for momentum, water, heat, and CO2
- Canopy boundary-layer conductance (Rb and kB-1 parameter; empirical and physically-based models)
- Monin-Obhukov length
- stability parameter zeta
- stability correction functions (different versions)
- roughness length for momentum (z0m) and heat (z0h)
- decoupling coefficient 'omega'
- Reynolds number
- wind speed at a given height from the logarithmic wind profile equation

## Surface conditions
- vapor pressure, specific humidity, and VPD at the big-leaf surface
- CO2 concentration at the big-leaf surface
- aerodynamic surface temperature
- radiometric surface temperature 

## Evapotranspiration (ET) and water-use efficiency (WUE)
- potential ET (Priestley-Taylor equation)
- reference ET (Penman-Monteith equation)
- imposed and equilibrum ET
- WUE, inherent WUE, underlying WUE

## Physiological variables
- canopy conductance (inverted Penman-Monteith equation)
- bulk intercellular CO2 concentration (Ci)
- bulk photosynthetic capacity (Vcmax25 and Jmax25)
- stomatal slope g1 (USO, Ball&Berry, Leuning models)
- stomatal sensitivity to VPD
- ecosystem light response, light-use efficiency

## Energy balance
- energy balance closure (EBR and slope method)
- biochemical energy
- energy-use efficiency

## Unit conversions
- conductance conversion from m s-1 to mol m-2 s-1 and vice versa
- conversions between humidity measures (vapor pressure, specific humidity, relative humidity, and VPD)
- conversion between latent heat flux (W m-2) and evapotranspiration (kg m-2 s-1)
- conversion between radiation in W m-2 and umol m-2 s-1
- carbon fluxes from umol m-2 s-1 to g m-2 day-1

# Contact
For questions, remarks, and suggestions please contact jknauer@bgc-jena.mpg.de

