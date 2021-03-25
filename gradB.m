function [gradB_x, gradB_y, gradB_z, gradBmod] = gradB(x,y,z)
% ------------------------------------------------------------------------
% gradB calculation
% ------------------------------------------------------------------------
  % Calculate  the gradient box
  dx = 0.01;	% grid step in X in m
  dy = 0.01;
  dz = 0.01;
  NB = 2;		% box order
  xx = linspace(-NB*dx, NB*dx, 2*NB+1);
  yy = linspace(-NB*dy, NB*dy, 2*NB+1);
  zz = linspace(-NB*dz, NB*dz, 2*NB+1);
  K = length(xx);
  M = length(yy);
  N = length(zz);
  [X,Y,Z] = meshgrid(xx,yy,zz);
    
  BX = zeros(K,M,N);
  BY = zeros(K,M,N);
  BZ = zeros(K,M,N);
  BB = zeros(K,M,N);
  [BX, BY, BZ, BB] = BfieldFast(X + x, Y + y, Z + z);

  % Calculate the gradient of  the module
  [dxB, dyB, dzB] = gradient(BB, dx, dy, dz); 
  gradB_x = dxB(NB+1, NB+1, NB+1);
  gradB_y = dyB(NB+1, NB+1, NB+1);
  gradB_z = dzB(NB+1, NB+1, NB+1);
  gradBmod = sqrt(gradB_x.^2 + gradB_y.^2 + gradB_z.^2);   
endfunction
