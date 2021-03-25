function GCdot = Lorentz_forceGC(GC)

% ------------------------------------------------------------------------
% LSODE equations for the Guiding Centre Approximation
% ------------------------------------------------------------------------
% Function that defines the equation of motion in
% terms of the 6 dimensional vector guding centre

  %Define the global variables
  global m q;
 
global VV = particle.guidingcenter.gc.position.v0;
global EFIT_R
global EFIT_Z
global EFIT_BPR 
global EFIT_BPZ
global EFIT_BT 

  % Initialize the solution array
  GCdot = zeros(4,1);
  
  % Calculate the B field and its gradient
  [Bx, By, Bz, B] = Bfield(GC(1),GC(2),GC(3));
  [gradB_x, gradB_y, gradB_z, gradBmod] = gradB(GC(1),GC(2),GC(3));
  
  % Set of equations for the velocity
  GCdot(1) = 0.5*((m*VV^2)/(q*B^3))*(1+(GC(4)/VV)^2)*(By*gradB_z - Bz*gradB_y) + GC(4)*Bx/B;
  GCdot(2) = 0.5*((m*VV^2)/(q*B^3))*(1+(GC(4)/VV)^2)*(Bz*gradB_x - Bx*gradB_z) + GC(4)*By/B;
  GCdot(3) = 0.5*((m*VV^2)/(q*B^3))*(1+(GC(4)/VV)^2)*(Bx*gradB_y - By*gradB_x) + GC(4)*Bz/B;
  
  % Set of equations for the acceleration
  GCdot(4) = -0.5*(VV^2 - GC(4)^2)*(Bx*gradB_x + By*gradB_y + Bz*gradB_z)/(B^2);
