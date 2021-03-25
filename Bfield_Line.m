function dBFLdot = Bfield_Line(P)
% ------------------------------------------------------------------------
% LSODE equations
% ------------------------------------------------------------------------
% Function that defines the equation of the field lines in terms of the
% particle position P

  %Define the global variables
  global EFIT_Rmax EFIT_Zmax EFIT_Zmin;
  
  % Check if the field line is outside EFIT boundary box:
  if (sqrt(P(1)^2 + P(2)^2) > EFIT_Rmax || P(3) > EFIT_Zmax || P(3) < EFIT_Zmin)
    return;
  endif

  % Initialize the solution array
  dBFLdot(1) = P(1);
  dBFLdot(2) = P(2);
  dBFLdot(3) = P(3);
  
  % Calculate the B field
  [Bx, By, Bz, BB] = BfieldFast(P(1),P(2),P(3));
  
  % Set of equations for the velocity
  dBFLdot(1) = Bx/BB;
  dBFLdot(2) = By/BB;
  dBFLdot(3) = Bz/BB;
  
endfunction
