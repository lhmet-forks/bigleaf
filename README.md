# README #

This repository serves as comprehensive collection of functions to calcultate additional variables from water, carbon and energy fluxes measured by EC data as well as meteorological data. The functions assume the vegetation to behave like a "big leaf", i.e. vertical variations within the canopy are ignored.


an (incomplete) list of what should be in the package:
# atmospheric variables:
- air density
- virtual temperature
- pressure from altitude
- psychrometric constant
- latent heat of vaporization
- slope of saturation vapor pressure curve
- saturation vapor pressure
- humidity conversions 

# physiological variables:
- Canopy conductance (PM, diffusion equation)
- Ci
- Gc conversion ms-1 to mol m-2 s-1

# ET
- PET, Priestley-Taylor

# aerodynamic properties:
- aerodynamic conductance
- quasi-laminar boundary layer conductance (kB)
- Monin-Obhukov length
- omega (different versions)
- Reynolds number

# surface conditions
- D0, T0, C0