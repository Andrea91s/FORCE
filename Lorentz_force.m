function PVdot = Lorentz_force(PV)
% ------------------------------------------------------------------------
% LSODE equations
% ------------------------------------------------------------------------
% Function that defines the equation of motion in
% terms of the 6 dimensional vector PV: Position and Velocity

  %Define the global variables
  global m q;
  global EFIT_Rmax EFIT_Zmax EFIT_Zmin;
  global EFIT_Raxis EFIT_BT0;
  
  % Initialize the solution array
  % PV = (x, y, z, vx, vy, vz)
  % PVdot = (dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt);
  PVdot = zeros(6,1);
  
  % Check if particle is outside EFIT boundary box:
  if (sqrt(PV(1)^2 + PV(2)^2) > EFIT_Rmax || PV(3) > EFIT_Zmax || PV(3) < EFIT_Zmin)
    return;
  endif
  
  
  % Calculate the B field
  [Bx, By, Bz] = BfieldFast(PV(1),PV(2),PV(3));
  
  if (EFIT_BT0 != 0)
  
    % Calculate the  TF ripple magnetic field
    [TFBx, TFBy, TFBz] = TBfieldRipple(PV(1), PV(2), PV(3), EFIT_Raxis, EFIT_BT0);
    
    % And the total magnetic field
    Bx = Bx + TFBx;
    By = By + TFBy;
    Bz = Bz + TFBz;
  
  endif
  
  % Calculate the E field
  [Ex, Ey, Ez] = Efield(PV(1),PV(2),PV(3));
  
  % Set of equations for the velocity
  PVdot(1) = PV(4);
  PVdot(2) = PV(5);
  PVdot(3) = PV(6);
  
  % Set of equations for the acceleration
  PVdot(4) = (PV(5)*Bz - PV(6)*By + Ex)*(q/m);
  PVdot(5) = (PV(6)*Bx - PV(4)*Bz + Ey)*(q/m);
  PVdot(6) = (PV(4)*By - PV(5)*Bx + Ez)*(q/m);
endfunction
