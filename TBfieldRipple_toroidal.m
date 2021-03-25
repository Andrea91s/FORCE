function [TFBR, TFBz, TFBphi, TFB] = TBfieldRipple_toroidal(R,Z,phi, Ro, Bto)
% Function that calculates the TF ripple magnetic field

% TF ripple parameters
Nc = 12;        % number of TF coils on MAST
Rc = 2.09;      % TF coils radius (m)

% Coordinates
% TF ripple magnetic field (cylindrical Coordinates):
TFBphi = Bto*(Ro./R).*((R/Rc).^Nc).*cos(Nc*phi);
TFBR = Bto*(Ro./R).*((R/Rc).^Nc).*sin(Nc*phi);
TFBz = zeros(size(R));

% Total Toroidal Field Ripple
TFB = sqrt(TFBR.^2 + TFBz.^2 + TFBphi.^2);

