function [TFBx, TFBy, TFBz, TFB, TFBphi, TFBR] = TBfieldRipple(x,y,z, Ro, Bto)
% Function that calculates the TF ripple magnetic field

% TF ripple parameters
Nc = 12;        % number of TF coils on MAST
Rc = 2.09;      % TF coils radius (m)

% Coordinates
R = sqrt(x.^2 + y.^2);
phi = atan2(y,x);

% TF ripple magnetic field (cylindircal Coordinates):
TFBphi = Bto*(Ro./R).*((R/Rc).^Nc).*cos(Nc*phi);
TFBR = Bto*(Ro./R).*((R/Rc).^Nc).*sin(Nc*phi);

% TF ripple magnetic field (carthesian coordinates):
TFBx = TFBR.*cos(phi) + TFBphi.*cos(phi+pi/2);
TFBy = TFBR.*sin(phi) + TFBphi.*sin(phi+pi/2);
TFBz = zeros(size(x));

% Total Toroidal Field Ripple
TFB = sqrt(TFBx.^2 + TFBy.^2 + TFBz.^2);

