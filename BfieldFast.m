% ------------------------------------------------------------------------
% B field in Tesla (array definition)
% ------------------------------------------------------------------------
function [BX, BY, BZ, BB] = BfieldFast(x,y,z)
global EFIT_R
global EFIT_Z
global EFIT_BPR 
global EFIT_BPZ
global EFIT_BT

	  %B = zeros(length(x),3);
	  r   = sqrt(x.^2 + y.^2);
	  bpr = interp2(EFIT_R, EFIT_Z, EFIT_BPR, r, z);
	  bpz = interp2(EFIT_R, EFIT_Z, EFIT_BPZ, r, z);
	  bt  = interp2(EFIT_R, EFIT_Z, EFIT_BT, r, z);
	  phi = atan2(y,x);
	  BX  = bpr.*cos(phi) + bt.*cos(phi+pi/2);
	  BY  = bpr.*sin(phi) + bt.*sin(phi+pi/2);
	  BZ  = bpz;	  
	  BB  = sqrt(BX.^2+BY.^2+BZ.^2);
endfunction	
