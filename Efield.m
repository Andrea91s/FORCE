function [Ex, Ey, Ez, Em] = Efield(x, y, z)
% Function that defines the E field in which the
% charged particle is moving.
% Units are in Volt/m 

% --------------------------------------------------	  
% Electric field
% --------------------------------------------------

	  Ex = zeros(size(x));
	  Ey = zeros(size(y));
	  Ez = zeros(size(z));
	  %Ez = 1E5*ones(size(z));
    Em = sqrt(Ex.^2+Ey.^2+Ez.^2);   
      
