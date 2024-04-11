Example parts:
1. The 'Very basic hydrologic' model (VBH): model description (Test model desciption.pdf), model MATLAB script (VBH_script.m)
		Run the script to generate example fluxes and stores
2. Linker example: data formatter function for the VBH model (data_repacker.m) and the other data used (Other data.txt)
		Run the repacker function on the example fluxes


MAITsim_18O should run on the outputs from the repacker function

*The VBH and the data repacker can be adapted to test the function of MAITsim 
(change flux magnitudes and storage in the VBH script; E/ET split, initial isotope compositions and atmosphere in the data repacker)