# README #

This repository serves as comprehensive collection of functions that can be used to calculate additional variables from EC flux data and accompanying meteorological measurements. All calculations are based on a 'big-leaf' representation of the the vegetation, i.e. vertical meteorological variations within the canopy are ignored and all fluxes are assumed to originate from a single plane at a certain height in the canopy.

It's supposed to develop into an R-package to be put on CRAN.


a list of functions that can be found in the package:

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


# unit conversions
- conductance conversion from ms-1 to mol m-2 s-1
- conversions between humidity measures vapor pressure, specific humidity, relative humidity, and VPD
- conversion between laten heat flux and evapotranspiration

# physiological variables:
- Canopy conductance (inverted PM)
- Bulk intercellulary CO2 concentrtation (Ci)


# ET and WUE
- Potential ET (Priestley-Taylor)
- imposed and equilibrum ET
- WUE, inherent WUE, underlying WUE
- g1 (USO, Ball&Berry, Leuning)


# aerodynamic properties:
- aerodynamic conductance (different versions) for water, heat, momentum, and CO2
- Canopy boundary-layer conductance (Rb and kB-1 parameter; different models)
- Monin-Obhukov length
- stability parameter zeta
- stability correction functions (different versions)
- roughness length for momentum (z0m) estimation
- omega (amphi- and hypostomatous formulation)
- Reynolds number
- wind speed at given height from wind profile equation


# surface conditions
- VPD, Ca, Temperature, e, q at the big-leaf surface

# energy balance
- biochemical energy
- energy balance ratio (EBR)
- "missing" energy