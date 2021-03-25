% Estimation of the NBI pitch angle on MAST
%
% This is done assuming that the NBI is injected with
% a tangency radius of 0.7 m parallel to the X axis this although
% not true does not matter since the magnetic field from EFIT is
% assumed to be toroidally symmetric.
%
% An additional hypthesys is that the initial velocity of the beam 
% ions is in the same direction as the injected neutral, that is
% parallel to the X axis. 

% Load an EFIT equilibria:
load EFIT_29881_at_0.250_s.mat 
MOrbit_global_variables;

% NBI tangency radius
RT = 0.7;

% Determines the coordinates of the beam path in the plasma
R = 2;
xm = sqrt(R^2 - RT^2);

% Defines the trajectory of the beam ions
x = linspace(-xm, xm, 1000);
y = RT;
z = 0;
r = sqrt(x.^2 + y^2);

% Evaluates the magnetic field along the NBI trajectory
[BX, BY, BZ, BB] = BfieldFast(x,y,z);

% NBI velocity versors
vx = ones(size(BX));
vy = zeros(size(BY));
vz = zeros(size(BZ));

% Evaluates the pitch angle
[mah, boh, pa] = projection([vx' vy' vz'], [BX' BY' BZ']);
clear mah boh 

% Evaluates lambda
lambda = cos(pa);

% Makes nice plots
figure(1, 'position', [100, 100, 600, 900])

subplot(3,1,1)
plot(r, BB, 'linewidth', 2)
axis([0.9*RT 2])
xlabel('Radial position (m)')
ylabel('Magnetic field |B| (Tesla)')
title(['B field along the NBI injection direction for R_T = ' num2str(RT) ' m'], 'fontsize', 12)

subplot(3,1,2)
plot(r, pa, 'linewidth', 2)
axis([0.9*RT 2])
xlabel('Radial position (m)')
ylabel('Pitch angle (rad)')
title(['Pitch-angle along the NBI injection direction for R_T = ' num2str(RT) ' m'], 'fontsize', 12)

subplot(3,1,3)
plot(r, lambda, 'linewidth', 2)
axis([0.9*RT 2])
xlabel('Radial position (m)')
ylabel('\lambda')
title(['Cosine of the pitch-angle along the NBI injection direction for R_T = ' num2str(RT) ' m'], 'fontsize', 12)
