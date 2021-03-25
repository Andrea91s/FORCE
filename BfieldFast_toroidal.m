% ------------------------------------------------------------------------
% B field in Tesla (array definition)
% ------------------------------------------------------------------------
function [bpr, bpz, bt, btok] = BfieldFast_toroidal(r,z,phi)
global EFIT_R
global EFIT_Z
global EFIT_BPR 
global EFIT_BPZ
global EFIT_BT
	  bpr = interp2(EFIT_R, EFIT_Z, EFIT_BPR, r, z);
	  bpz = interp2(EFIT_R, EFIT_Z, EFIT_BPZ, r, z);
	  bt = interp2(EFIT_R, EFIT_Z, EFIT_BT, r, z);	  
	  btok = sqrt(bpr.^2 + bpz.^2 + bt.^2);
endfunction	
