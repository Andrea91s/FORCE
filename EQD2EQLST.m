function [EQLST] = EQD2EQLST(EQD)
% Function that convert the structure EQD from EQDSK files to the MORBIT EQLST


EQLST.R = EQD.R_grid;
EQLST.Z = EQD.Z_grid;
EQLST.magneticfield.polodailBr = EQD.B.POL.R;
EQLST.magneticfield.polodailBz = EQD.B.POL.Z;
EQLST.magneticfield.toroidalB = EQD.B.TOR;
EQLST.poloidalFlux = EQD.psi_2D;
EQLST.magneticAxisr = EQD.R_axis;
EQLST.magneticAxisz = EQD.Z_axis;
EQLST.poloidalFluxBoundary = EQD.psi_edge;
EQLST.rb = EQD.LCFS_R;
EQLST.zb = EQD.LCFS_Z;
EQLST.pulseNumber = EQD.pulseNumber;
EQLST.selectedtime = EQD.selectedtime;
