% ------------------------------------------------------------------------
% B field in Tesla
% ------------------------------------------------------------------------
function [B] = Bfield(x,y,z)
global EFIT_R
global EFIT_Z
global EFIT_BPR 
global EFIT_BPZ
global EFIT_BT

	  %B = zeros(1,3);
	  r   = sqrt(x^2 + y^2);
	  bpr = interp2(EFIT_R, EFIT_Z, EFIT_BPR, r, z);
	  bpz = interp2(EFIT_R, EFIT_Z, EFIT_BPZ, r, z);
	  bt  = interp2(EFIT_R, EFIT_Z, EFIT_BT, r, z);
	  phi = atan2(y,x);
	  Bx  = bpr*cos(phi) + bt*cos(phi+pi/2);
	  By  = bpr*sin(phi) + bt*sin(phi+pi/2);
	  Bz  = bpz;

	  B(1,1) = Bx;
	  B(1,2) = By;
	  B(1,3) = Bz;
	  
endfunction	
