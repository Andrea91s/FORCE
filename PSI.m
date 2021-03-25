% Poloidal Flux
function [PSI] = PSI(R,z)

global EFIT_R
global EFIT_Z
global EFIT_PSI

	  %PSI = interp2(EFIT_R, EFIT_Z, EFIT_PSI, R, z);
	  PSI = interp2(EFIT_R, EFIT_Z, EFIT_PSI, R, z, 'pchip');
  
endfunction	
