% Definition of global variables required by different scripts according to the table below:
%
% variable  Bfield.m    Lorentz_force.m     BfieldFast.m    MOrbit_post_processor   PSI.m
% ---------------------------------------------------------------------------------------
% EFIT_R    X                               X                                       X
% EFIT_Z    X                               X                                       X
% EFIT_BPR  X                               X
% EFIT_BPZ  X                               X
% EFIT_BT   X                               X
% EFIT_PSI                                                                          X
% EFIT_Rmin             X
% EFIT_Rmax             X
% EFIT_Zmin             X
% EFIT_Zmax             X
% EFIT_BB0                                                  X
% EFIT_PSIa                                                 X
% EFIT_BT0              X (set to 0 to run without TF ripples)

global mu0 = 4*pi*1E-7;
global EFIT_R = EQLST.R;
global EFIT_Z = EQLST.Z;
global EFIT_BPR = EQLST.magneticfield.polodailBr;
global EFIT_BPZ = EQLST.magneticfield.polodailBz;
global EFIT_BT = EQLST.magneticfield.toroidalB;
global EFIT_PSI = EQLST.poloidalFlux;
global EFIT_Rmin = min(EQLST.R(:));
global EFIT_Rmax =  max(EQLST.R(:));
global EFIT_Zmin =  min(EQLST.Z(:));
global EFIT_Zmax = max(EQLST.Z(:));
global EFIT_Raxis = EQLST.magneticAxisr;

% Evaluates the magnetic field on axis and the edge
[BX, BY, BZ, BB0] = BfieldFast(EQLST.magneticAxisr, 0, EQLST.magneticAxisz);
global EFIT_BB0 = BB0;
%global EFIT_BT0 = BY;
global EFIT_BT0 = 0;
global EFIT_PSIa = EQLST.poloidalFluxBoundary;
clear BX BY BZ BB0
        
        
