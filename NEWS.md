NEWS / ChangeLog for bigleaf
---------------------------------

0.4.4    13.03.2018   
         - reference.ET: argument Gs renamed to Gs_ref
		  
0.4.3    12.03.2018
         - bigleaf.constants: restructured in a way that allows changing the constants for each function call
		  
0.4.2    12.03.2018
         - Esat.slope: constants (a,b,c) from Allen et al. 1998 added as option. Argument "Esat.formula" added to every 
		               function that calculates esat or slope of esat to allow consistency among functions.  

		 - light.response: argument "..." added to the nls function
         - stomatal.slope: argument "..." within each call of nls in the function	

0.4.1    08.03.2018
         - aerodynamic.conductance: Ga_CO2 added to function output 		 
		 
0.4.0    08.03.2018
         - aerodynamic.conductance: Rb can be calculated for other quantities if the respective Schmidt number is provided.
                                    I.e. new arguments "Sc" and "Sc_name" as the value of the Schmidt number, and the name
                                    of the quantity for which Sc is provided, respectively.
         - bigleaf.constants: Prandtl number (Pr) and Schmidt number for CO2 (Sc_CO2) added
		 - light.response: bug fix: +Reco replaced by -Reco (sign was reversed)
		 
0.3.2    06.03.2018	
         - rH.to.VPD: if statement vectorized

0.3.2    12.02.2018
         - radiometric.surface.temp: Output renamed	("." replaced by "_")	

0.3.1    09.02.2018
         - Monin.Obukhov.length: default arguments added
		 
         		 