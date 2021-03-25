function [particle] = MOrbit_post_processor(particle)
% ------------------------------------------------------------------------
% FULL ORBIT POST PROCESSOR
% ------------------------------------------------------------------------
tic();
printf('%s', 'Full orbit post-processor: '); fflush(stdout());

% Declares global variables
global EFIT_BB0;
global EFIT_PSIa;
global EFIT_Raxis EFIT_BT0;


% Compute solutions
particle.orbit.velocity.v = sqrt(particle.orbit.velocity.vx.^2 + particle.orbit.velocity.vy.^2 + particle.orbit.velocity.vz.^2);
particle.orbit.velocity.vxy = sqrt(particle.orbit.velocity.vx.^2 + particle.orbit.velocity.vy.^2);

% Calculates the maximum distance travelled in a single step
particle.orbit.max_distance = max(particle.orbit.velocity.v)*particle.times.dt;

% Radial coordinate (for poloidal projection of the orbits)
particle.orbit.position.R = sqrt(particle.orbit.position.x.^2 + particle.orbit.position.y.^2);

% Toroidal angle
particle.orbit.position.phi = atan2(particle.orbit.position.y, particle.orbit.position.x);
particle.orbit.position.phi = unwrap(particle.orbit.position.phi);

% Calculates the Toroidal Velocity
dt = particle.orbit.t(2) - particle.orbit.t(1);
particle.orbit.velocity.vtor = particle.orbit.position.R.*gradient(particle.orbit.position.phi, dt);

%Compute the elextric field along the trajectory
[particle.orbit.electric_field.Ex, particle.orbit.electric_field.Ey, particle.orbit.electric_field.Ez] = Efield(particle.orbit.position.x, particle.orbit.position.y, particle.orbit.position.z);
particle.orbit.electric_field.E = sqrt(particle.orbit.electric_field.Ex.^2 + particle.orbit.electric_field.Ey.^2 + particle.orbit.electric_field.Ez.^2);


% Compute the electric potential alonng the trajectory
% using method 1 in Apostol Calculus III page 298 and 297
[EX0, dummy, dummy] = Efield(particle.orbit.position.x, particle.orbit.position.y(1), particle.orbit.position.z(1));
[dummy, EY0, dummy] = Efield(particle.orbit.position.x(1), particle.orbit.position.y, particle.orbit.position.z(1));
[dummy, dummy, EZ0] = Efield(particle.orbit.position.x(1), particle.orbit.position.y(1), particle.orbit.position.z);
F1 = -cumtrapz(particle.orbit.position.x, EX0);
F2 = -cumtrapz(particle.orbit.position.y, EY0);
F3 = -cumtrapz(particle.orbit.position.z, EZ0);
particle.orbit.electric_field.phi = F1 + F2 + F3;


% Calculate the particle energy
particle.orbit.energy.kinetic = 0.5*particle.mass*particle.orbit.velocity.v.^2;
particle.orbit.energy.potential = particle.charge*particle.orbit.electric_field.phi;
particle.orbit.energy.total = particle.orbit.energy.kinetic + particle.orbit.energy.potential;

% Compute B fields along the trajectory
[AX, AY, AZ, AA] = BfieldFast(particle.orbit.position.x, particle.orbit.position.y, particle.orbit.position.z); 
[TFBx, TFBy, TFBz] = TBfieldRipple(particle.orbit.position.x, particle.orbit.position.y, particle.orbit.position.z, EFIT_Raxis, EFIT_BT0); 
AX = AX + TFBx;
AY = AY + TFBy;
AZ = AZ + TFBz;
AA = sqrt(AX.^2 + AY.^2 + AZ.^2);
particle.orbit.magnetic_field.Bx = AX;
particle.orbit.magnetic_field.By = AY;
particle.orbit.magnetic_field.Bz = AZ;
particle.orbit.magnetic_field.B = AA;

% Compute B fields along the R,Z and phi
[TR, TZ, TPHI, TT] = BfieldFast_toroidal(particle.orbit.position.R, particle.orbit.position.z, particle.orbit.position.phi); 
[TFBR, TFBZ, TFBPHI] = TBfieldRipple_toroidal(particle.orbit.position.R, particle.orbit.position.z, particle.orbit.position.phi, EFIT_Raxis, EFIT_BT0); 
TR = TR + TFBR;
TZ = TZ + TFBZ;
TPHI = TPHI + TFBPHI;
TT = sqrt(TR.^2 + TZ.^2 + TPHI.^2);
particle.orbit.magnetic_field.BR = TR;
particle.orbit.magnetic_field.BZ = TZ;
particle.orbit.magnetic_field.BPHI = TPHI;
particle.orbit.magnetic_field.BT = TT;
clear TR TZ TPHI TT TFBR TFBZ TFBPHI

% Calculate the magnetic field versors
particle.orbit.magnetic_field.bx = particle.orbit.magnetic_field.Bx./particle.orbit.magnetic_field.B;
particle.orbit.magnetic_field.by = particle.orbit.magnetic_field.By./particle.orbit.magnetic_field.B;
particle.orbit.magnetic_field.bz = particle.orbit.magnetic_field.Bz./particle.orbit.magnetic_field.B;
clear AX AY AZ AA TFBx TFBy TFBz


% Compute v parallel
vdotB = particle.orbit.velocity.vx.*particle.orbit.magnetic_field.bx + particle.orbit.velocity.vy.*particle.orbit.magnetic_field.by + particle.orbit.velocity.vz.*particle.orbit.magnetic_field.bz;
particle.orbit.velocity.v_pll_x = vdotB.*particle.orbit.magnetic_field.bx;
particle.orbit.velocity.v_pll_y = vdotB.*particle.orbit.magnetic_field.by;
particle.orbit.velocity.v_pll_z = vdotB.*particle.orbit.magnetic_field.bz;
particle.orbit.velocity.v_pll = sqrt(particle.orbit.velocity.v_pll_x.^2 + particle.orbit.velocity.v_pll_y.^2 + particle.orbit.velocity.v_pll_z.^2);
particle.orbit.energy.kinetic_parl = 0.5*particle.mass*particle.orbit.velocity.v_pll.^2;

% Compute v perp
particle.orbit.velocity.v_perp_x = particle.orbit.velocity.vx - particle.orbit.velocity.v_pll_x;
particle.orbit.velocity.v_perp_y = particle.orbit.velocity.vy - particle.orbit.velocity.v_pll_y;
particle.orbit.velocity.v_perp_z = particle.orbit.velocity.vz - particle.orbit.velocity.v_pll_z;
particle.orbit.velocity.v_perp = sqrt(particle.orbit.velocity.v_perp_x.^2 + particle.orbit.velocity.v_perp_y.^2 + particle.orbit.velocity.v_perp_z.^2);
particle.orbit.energy.kinetic_perp = 0.5*particle.mass*particle.orbit.velocity.v_perp.^2;


% Compute pitch angle
%[dummy1, dummy2, dummy3] = projection(V1,V2);
%particle.orbit.pitch_angle = dummy3;
%clear dummy1 dummy2 dummy3 V1 V2
V1 = [particle.orbit.velocity.vx, particle.orbit.velocity.vy, particle.orbit.velocity.vz];
V2 = [particle.orbit.magnetic_field.Bx, particle.orbit.magnetic_field.By, particle.orbit.magnetic_field.Bz];
particle.orbit.pitch_angle = vectorAngle3d(V1, V2);
clear V1 V2

% Computes Lambda
particle.orbit.lambda = cos(particle.orbit.pitch_angle);

% Compute the rotation angle around B
particle.gyro_orbit.theta = atan2(particle.orbit.velocity.v_perp_y, particle.orbit.velocity.v_perp_x);
particle.gyro_orbit.theta_unwrapped = unwrap(particle.gyro_orbit.theta);


% Larmor radius
% radius of gyration
% gyro-radius
% cyclotron radius
particle.orbit.larmor_radius = (particle.mass/abs(particle.charge))*particle.orbit.velocity.v_perp./particle.orbit.magnetic_field.B;

% Angular velocity magnitude is known: 
% Angular frequency of gyration
% AKA ciclotron frequency
% AKA Larmor frequency
% AKA gyro-frequency
particle.orbit.angular_velocity = particle.charge*particle.orbit.magnetic_field.B/particle.mass;
particle.orbit.angular_frequency = particle.orbit.angular_velocity/(2*pi);

% magnetic moment
particle.orbit.magnetic_moment.mu = particle.mass*(particle.orbit.velocity.v_perp.^2)./(2*particle.orbit.magnetic_field.B);
particle.orbit.magnetic_moment.mux = -particle.orbit.magnetic_moment.mu.*particle.orbit.magnetic_field.bx;
particle.orbit.magnetic_moment.muy = -particle.orbit.magnetic_moment.mu.*particle.orbit.magnetic_field.by;
particle.orbit.magnetic_moment.muz = -particle.orbit.magnetic_moment.mu.*particle.orbit.magnetic_field.bz;


% Calculates the toroidal versor
%particle.orbit.position.tor_versor.x = cos(pi/2).*particle.orbit.position.x - sin(pi/2).*particle.orbit.position.y;
%particle.orbit.position.tor_versor.y = sin(pi/2).*particle.orbit.position.x + cos(pi/2).*particle.orbit.position.y;
%particle.orbit.position.tor_versor.z = zeros(size(particle.orbit.position.tor_versor.x));
%particle.orbit.position.tor_versor.norm = sqrt(particle.orbit.position.tor_versor.x.^2 + particle.orbit.position.tor_versor.y.^2 + particle.orbit.position.tor_versor.z.^2);
%particle.orbit.position.tor_versor.ux = particle.orbit.position.tor_versor.x./particle.orbit.position.tor_versor.norm;
%particle.orbit.position.tor_versor.uy = particle.orbit.position.tor_versor.y./particle.orbit.position.tor_versor.norm;
%particle.orbit.position.tor_versor.uz = particle.orbit.position.tor_versor.z./particle.orbit.position.tor_versor.norm;
%particle.orbit.position.tor_versor.vtx = (particle.orbit.velocity.vx.*particle.orbit.position.tor_versor.ux + ...
%                                          particle.orbit.velocity.vy.*particle.orbit.position.tor_versor.uy + ...
%                                          particle.orbit.velocity.vz.*particle.orbit.position.tor_versor.uz).*particle.orbit.position.tor_versor.ux;
%particle.orbit.position.tor_versor.vty = (particle.orbit.velocity.vx.*particle.orbit.position.tor_versor.ux + ...
%                                          particle.orbit.velocity.vy.*particle.orbit.position.tor_versor.uy + ...
%                                          particle.orbit.velocity.vz.*particle.orbit.position.tor_versor.uz).*particle.orbit.position.tor_versor.uy;
%particle.orbit.position.tor_versor.vtz = (particle.orbit.velocity.vx.*particle.orbit.position.tor_versor.ux + ...
%                                          particle.orbit.velocity.vy.*particle.orbit.position.tor_versor.uy + ...
%                                          particle.orbit.velocity.vz.*particle.orbit.position.tor_versor.uz).*particle.orbit.position.tor_versor.uz;  
%particle.orbit.position.tor_versor.vtor = sqrt(particle.orbit.position.tor_versor.vtx.^2 + particle.orbit.position.tor_versor.vty.^2 + particle.orbit.position.tor_versor.vtz.^2);                                          

%particle.orbit.velocity.vradial =  (particle.orbit.velocity.vx.*particle.orbit.position.x + particle.orbit.velocity.vy.*particle.orbit.position.y)./(particle.orbit.position.R.^2);
%particle.orbit.velocity.vtoroidal = particle.orbit.velocity.vx.*particle.orbit.position.tor_versor.ux + particle.orbit.velocity.vy.*particle.orbit.position.tor_versor.uy;

% XY components of the radial velocity
particle.orbit.velocity.vrx = (particle.orbit.velocity.vx.*particle.orbit.position.x + particle.orbit.velocity.vy.*particle.orbit.position.y).*particle.orbit.position.x./(particle.orbit.position.R.^2);
particle.orbit.velocity.vry = (particle.orbit.velocity.vx.*particle.orbit.position.x + particle.orbit.velocity.vy.*particle.orbit.position.y).*particle.orbit.position.y./(particle.orbit.position.R.^2);
particle.orbit.velocity.vr = sqrt(particle.orbit.velocity.vrx.^2 + particle.orbit.velocity.vry.^2);

% XY components of the toroidal velocity
particle.orbit.velocity.vtx = (particle.orbit.velocity.vy.*particle.orbit.position.x - particle.orbit.velocity.vx.*particle.orbit.position.y).*(-particle.orbit.position.y)./(particle.orbit.position.R.^2);
particle.orbit.velocity.vty = (particle.orbit.velocity.vy.*particle.orbit.position.x - particle.orbit.velocity.vx.*particle.orbit.position.y).*particle.orbit.position.x./(particle.orbit.position.R.^2);
particle.orbit.velocity.vt = sqrt(particle.orbit.velocity.vtx.^2 + particle.orbit.velocity.vty.^2);

% Module of the velocity in the XY plane from its toroidal and radial components which is equal to sqrt(vx^2+vy^2);
particle.orbit.velocity.vrt = sqrt(particle.orbit.velocity.vr.^2 + particle.orbit.velocity.vt.^2);

% Calculation of the components of vxy in the radial and toroidal directions
V1 = [particle.orbit.velocity.vx, particle.orbit.velocity.vy];
V2 = [particle.orbit.velocity.vrx, particle.orbit.velocity.vry];
particle.orbit.velocity.angle = vectorAngle(V1, V2);
clear V1 V2
particle.orbit.velocity.vxyrad = particle.orbit.velocity.vxy.*cos(particle.orbit.velocity.angle);
particle.orbit.velocity.vxytor = particle.orbit.velocity.vxy.*sin(particle.orbit.velocity.angle);

% Calculation again of the total xy velocity from the above components for double check
particle.orbit.velocity.vxyver = sqrt(particle.orbit.velocity.vxyrad.^2 + particle.orbit.velocity.vxytor.^2);

% And an even simpler toroidal velocity definition
particle.orbit.velocity.vphi = -particle.orbit.velocity.vx.*sin(particle.orbit.position.phi) + particle.orbit.velocity.vy.*cos(particle.orbit.position.phi);


% Toroidal Canonical Angular Momentum
particle.orbit.PSI = PSI(particle.orbit.position.R, particle.orbit.position.z);
%particle.orbit.Ptor = particle.mass*particle.orbit.position.R.*particle.orbit.velocity.vtor - particle.charge*particle.orbit.PSI;
particle.orbit.Ptor = particle.mass*particle.orbit.position.R.*particle.orbit.velocity.vphi + particle.charge*particle.orbit.PSI;

% Calculate LAMBDA and the normalized toroidal canonical momentum
particle.orbit.LAMBDA = particle.orbit.magnetic_moment.mu*EFIT_BB0./(particle.orbit.energy.total);
particle.orbit.PtorPSIa = (particle.orbit.Ptor/physical_constant('electron volt'))/EFIT_PSIa;

printf('%f s\n', toc()); fflush(stdout());

