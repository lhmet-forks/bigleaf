# README #

This repository serves as comprehensive collection of functions that can be used to calculate additional variables from EC flux data and accompanying meteorological data. The functions assume that the vegetation behaves like a "big leaf", i.e. vertical variations within the canopy are ignored.


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
- aerodynamic conductance (different versions)
- quasi-laminar boundary layer conductance (kB)
- Monin-Obhukov length
- stability parameter zeta
- omega (different versions)
- Reynolds number
- stability correction functions (different versions)

# surface conditions
- D0, T0, C0

# energy balance
- Bowen ratio
- energy balance ratio (EBR)
- "missing" energy
- Gc uncertainty based on EBR