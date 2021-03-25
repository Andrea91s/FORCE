function [particle] = MOrbit_guiding_centre(particle)
% ------------------------------------------------------------------------
% GUIDING CENTRE CALCULATIONS from ORBIT calculations
% ------------------------------------------------------------------------

% Decalare global variables
global EFIT_Raxis EFIT_BT0;

tic();
printf('%s', 'Guiding centre calculations: '); fflush(stdout());
% Calculates the gyration radius vector:
%    calculation of v X B
particle.orbit.v_x_B_x = particle.orbit.velocity.vy.*particle.orbit.magnetic_field.Bz - particle.orbit.velocity.vz.*particle.orbit.magnetic_field.By;
particle.orbit.v_x_B_y = particle.orbit.velocity.vz.*particle.orbit.magnetic_field.Bx - particle.orbit.velocity.vx.*particle.orbit.magnetic_field.Bz;
particle.orbit.v_x_B_z = particle.orbit.velocity.vx.*particle.orbit.magnetic_field.By - particle.orbit.velocity.vy.*particle.orbit.magnetic_field.Bx;

%    calculation of |v X B|
particle.orbit.v_x_B_mod = sqrt(particle.orbit.v_x_B_x.^2 + particle.orbit.v_x_B_y.^2 + particle.orbit.v_x_B_z.^2);

%    calculation of  the gyration radius vector rg =  (m*Vperp/qB) (v X B)/|v x B|
particle.orbit.gyro_radius_x = (particle.mass*particle.orbit.velocity.v_perp./(particle.charge*particle.orbit.magnetic_field.B)).*(particle.orbit.v_x_B_x./particle.orbit.v_x_B_mod);
particle.orbit.gyro_radius_y = (particle.mass*particle.orbit.velocity.v_perp./(particle.charge*particle.orbit.magnetic_field.B)).*(particle.orbit.v_x_B_y./particle.orbit.v_x_B_mod);
particle.orbit.gyro_radius_z = (particle.mass*particle.orbit.velocity.v_perp./(particle.charge*particle.orbit.magnetic_field.B)).*(particle.orbit.v_x_B_z./particle.orbit.v_x_B_mod);
particle.orbit.gyro_radius_rmod =  sqrt(particle.orbit.gyro_radius_x.^2 + particle.orbit.gyro_radius_y.^2 + particle.orbit.gyro_radius_z.^2);

               
%    and then the guiding centre coordinates in the reference frame
particle.orbit.guidingcenter_x = particle.orbit.position.x + particle.orbit.gyro_radius_x;
particle.orbit.guidingcenter_y = particle.orbit.position.y + particle.orbit.gyro_radius_y;
particle.orbit.guidingcenter_z = particle.orbit.position.z + particle.orbit.gyro_radius_z;


% Radial coordinate (for poloidal projection of the orbits)
particle.orbit.guidingcenter_R = sqrt(particle.orbit.guidingcenter_x.^2 + particle.orbit.guidingcenter_y.^2);

% Toroidal angle
particle.orbit.guidingcenter_phi = atan2(particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_x);
particle.orbit.guidingcenter_phi = unwrap(particle.orbit.guidingcenter_phi);

% Calculates the Toroidal Velocity
dt = particle.orbit.t(2) - particle.orbit.t(1);
particle.orbit.guidingcenter_vtor = particle.orbit.guidingcenter_R.*gradient(particle.orbit.guidingcenter_phi, dt);

%   and its velocity
particle.orbit.guidingcenter_vx = gradient(particle.orbit.guidingcenter_x, dt);
particle.orbit.guidingcenter_vy = gradient(particle.orbit.guidingcenter_y, dt);
particle.orbit.guidingcenter_vz = gradient(particle.orbit.guidingcenter_z, dt);
particle.orbit.guidingcenter_v = particle.orbit.velocity.v;
%particle.orbit.guidingcenter_v = sqrt(particle.orbit.guidingcenter_vx.^2 + particle.orbit.guidingcenter_vy.^2 + particle.orbit.guidingcenter_vz.^2);
particle.orbit.guidingcenter_VX = particle.orbit.velocity.vx + (particle.mass/particle.charge)* gradient(particle.orbit.v_x_B_x./particle.orbit.magnetic_field.B.^2, dt);
particle.orbit.guidingcenter_VY = particle.orbit.velocity.vy + (particle.mass/particle.charge)* gradient(particle.orbit.v_x_B_y./particle.orbit.magnetic_field.B.^2, dt);
particle.orbit.guidingcenter_VZ = particle.orbit.velocity.vz + (particle.mass/particle.charge)* gradient(particle.orbit.v_x_B_z./particle.orbit.magnetic_field.B.^2, dt);
particle.orbit.guidingcenter_V = sqrt(particle.orbit.guidingcenter_VX.^2 + particle.orbit.guidingcenter_VY.^2 + particle.orbit.guidingcenter_VZ.^2);

% Calculating the B field at the guiding centre 
[particle.orbit.guidingcenter_Bx, particle.orbit.guidingcenter_By, particle.orbit.guidingcenter_Bz, particle.orbit.guidingcenter_B] = BfieldFast(particle.orbit.guidingcenter_x , particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_z);
[TFBx, TFBy, TFBz] = TBfieldRipple(particle.orbit.guidingcenter_x , particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_z, EFIT_Raxis, EFIT_BT0);
particle.orbit.guidingcenter_Bx = particle.orbit.guidingcenter_Bx + TFBx;
particle.orbit.guidingcenter_By = particle.orbit.guidingcenter_By + TFBy;
particle.orbit.guidingcenter_Bz = particle.orbit.guidingcenter_Bz + TFBz;
particle.orbit.guidingcenter_B = sqrt(particle.orbit.guidingcenter_Bx.^2 + particle.orbit.guidingcenter_By.^2 + particle.orbit.guidingcenter_Bz.^2);
clear TFBx TFBy TFBz
particle.orbit.guidingcenter_bx = particle.orbit.guidingcenter_Bx./ particle.orbit.guidingcenter_B;
particle.orbit.guidingcenter_by = particle.orbit.guidingcenter_By./ particle.orbit.guidingcenter_B;
particle.orbit.guidingcenter_bz = particle.orbit.guidingcenter_Bz./ particle.orbit.guidingcenter_B;

% Now I can calculate the difference (gradient) of the magnetic fields at the particle
% position and at the guiding centre position
particle.orbit.guidingcenter_deltaBx = particle.orbit.magnetic_field.Bx - particle.orbit.guidingcenter_Bx;
particle.orbit.guidingcenter_deltaBy = particle.orbit.magnetic_field.By - particle.orbit.guidingcenter_By;
particle.orbit.guidingcenter_deltaBz = particle.orbit.magnetic_field.Bz - particle.orbit.guidingcenter_Bz;
particle.orbit.guidingcenter_deltaB = sqrt(particle.orbit.guidingcenter_deltaBx.^2 + particle.orbit.guidingcenter_deltaBy.^2 + particle.orbit.guidingcenter_deltaBz.^2);
particle.orbit.guidingcenter_deltaBoverB = particle.orbit.guidingcenter_deltaB./particle.orbit.magnetic_field.B;



% Calculating v_par and v_per for the guiding centre
particle.orbit.guidingcenter.v_pll_x = (particle.orbit.guidingcenter_vx.*particle.orbit.guidingcenter_bx + ...
				      particle.orbit.guidingcenter_vy.*particle.orbit.guidingcenter_by + ...
				      particle.orbit.guidingcenter_vz.*particle.orbit.guidingcenter_bz).*particle.orbit.guidingcenter_bx;
particle.orbit.guidingcenter.v_pll_y = (particle.orbit.guidingcenter_vx.*particle.orbit.guidingcenter_bx + ...
				  particle.orbit.guidingcenter_vy.*particle.orbit.guidingcenter_by + ...
				  particle.orbit.guidingcenter_vz.*particle.orbit.guidingcenter_bz).*particle.orbit.guidingcenter_by;
particle.orbit.guidingcenter.v_pll_z = (particle.orbit.guidingcenter_vx.*particle.orbit.guidingcenter_bx + ...
				  particle.orbit.guidingcenter_vy.*particle.orbit.guidingcenter_by + ...
				  particle.orbit.guidingcenter_vz.*particle.orbit.guidingcenter_bz).*particle.orbit.guidingcenter_bz;
particle.orbit.guidingcenter.v_parl = sqrt(particle.orbit.guidingcenter.v_pll_x.^2 + particle.orbit.guidingcenter.v_pll_y.^2 + particle.orbit.guidingcenter.v_pll_z.^2);    
particle.orbit.guidingcenter.v_perp = sqrt(particle.orbit.guidingcenter_v.^2 - particle.orbit.guidingcenter.v_parl.^2);

% Calculating pitch angle and lambda for the guiding centre
V1 = [particle.orbit.guidingcenter_vx, particle.orbit.guidingcenter_vy, particle.orbit.guidingcenter_vx];
V2 = [particle.orbit.guidingcenter_Bx, particle.orbit.guidingcenter_By, particle.orbit.guidingcenter_Bz];
particle.orbit.guidingcenter.pitch_angle = vectorAngle3d(V1, V2);
clear V1 V2
particle.orbit.guidingcenter.lambda = cos(particle.orbit.guidingcenter.pitch_angle);

% Evaluates the parallel component of the velocity
particle.orbit.guidingcenter.v_parallel = particle.orbit.guidingcenter_v.*particle.orbit.guidingcenter.lambda;

% Evaluate the parallel kinetic energy
particle.orbit.guidingcenter.E_parallel = 0.5*particle.mass*particle.orbit.guidingcenter.v_parallel.^2;  
particle.orbit.guidingcenter.E_parallel = energy_unit (particle.orbit.guidingcenter.E_parallel, conversion = 0);


% Evaluates the magnetic moment at the centre og gyration                                   
particle.orbit.guidingcenter_mu = particle.mass*(particle.orbit.velocity.v_perp.^2)./(2*particle.orbit.guidingcenter_B);
particle.orbit.guidingcenter_rl = (particle.mass/abs(particle.charge))*particle.orbit.velocity.v_perp./particle.orbit.guidingcenter_B;
particle.orbit.guidingcenter_omega = particle.charge*particle.orbit.guidingcenter_B/particle.mass;
particle.orbit.guidingcenter_frequency = particle.orbit.guidingcenter_omega/(2*pi);

printf('%f s\n', toc()); fflush(stdout());


