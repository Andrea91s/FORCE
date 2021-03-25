function [particle] = MOrbit_guiding_center_approx(particle)
% ------------------------------------------------------------------------
% Initial conditions
% ------------------------------------------------------------------------
tic();
printf('%s', 'Guiding centre approximation solver: ');


% The initial coordinates of the guiding centre are obtained from the 
% full orbit motion

particle.guidingcenter.gc.position.x0 = particle.orbit.guidingcenter_x(1);
particle.guidingcenter.gc.position.y0 = particle.orbit.guidingcenter_y(1);
particle.guidingcenter.gc.position.z0 = particle.orbit.guidingcenter_z(1);
particle.guidingcenter.gc.position.v0 = particle.orbit.guidingcenter_v(1);
particle.guidingcenter.gc.position.v_par0 = particle.orbit.guidingcenter_v_pll(1);


% Calculates the magnetic field at the initial guiding centre position
[Bx0, By0, Bz0, B0] = BfieldFast(particle.guidingcenter.gc.position.x0, particle.guidingcenter.gc.position.y0, particle.guidingcenter.gc.position.z0);
particle.guidingcenter.gc.Bx0 = Bx0;
particle.guidingcenter.gc.By0 = By0;
particle.guidingcenter.gc.Bz0 = Bz0;
particle.guidingcenter.gc.B0 = B0;
clear Bx0 By0 Bz0 B0


% ------------------------------------------------------------------------
% Solve the GV euqation of motion
% ------------------------------------------------------------------------
[particle] = MOrbit_solveGC(particle);




% ------------------------------------------------------------------------
% Post processor
% ------------------------------------------------------------------------
particle.guidingcenter.gc.velocity.v_par = abs(particle.guidingcenter.gc.velocity.v_par);

% Radial coordinate (for poloidal projection of the orbits)
particle.guidingcenter.gc.position.R = sqrt(particle.guidingcenter.gc.position.x.^2 + particle.guidingcenter.gc.position.y.^2);

%   and its velocity
dtgc = particle.guidingcenter.gc.t(2) - particle.guidingcenter.gc.t(1);
particle.guidingcenter.gc.velocity.vx = gradient(particle.guidingcenter.gc.position.x, dtgc);
particle.guidingcenter.gc.velocity.vy = gradient(particle.guidingcenter.gc.position.y, dtgc);
particle.guidingcenter.gc.velocity.vz = gradient(particle.guidingcenter.gc.position.z, dtgc);
particle.guidingcenter.gc.velocity.vmod = sqrt(particle.guidingcenter.gc.velocity.vx.^2 + particle.guidingcenter.gc.velocity.vy.^2 + particle.guidingcenter.gc.velocity.vz.^2);


% Evaluates the B field at the guiding centre position
[particle.guidingcenter.gc.Bx, particle.guidingcenter.gc.By, particle.guidingcenter.gc.Bz, particle.guidingcenter.gc.B]  = Bfield(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, particle.guidingcenter.gc.position.z);
%particle.guidingcenter.gc.B = particle.guidingcenter.gc.B';


% Evaluates the E field at the guiding centre position
[particle.guidingcenter.gc.Ex, particle.guidingcenter.gc.Ey, particle.guidingcenter.gc.Ez, particle.guidingcenter.gc.E]  = Efield(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, particle.guidingcenter.gc.position.z);


% Compute the particle properties from the guiding centre
% position, such as v, E, mu ...
particle.guidingcenter.gc.velocity.v = particle.initialconditions.v0;
particle.guidingcenter.gc.velocity.v_per = sqrt(particle.guidingcenter.gc.velocity.v.^2 - particle.guidingcenter.gc.velocity.v_par.^2);
particle.guidingcenter.gc.magnetic_moment.mu = particle.mass*(particle.guidingcenter.gc.velocity.v_per.^2)./(2*particle.guidingcenter.gc.B);
particle.guidingcenter.gc.lambda = particle.guidingcenter.gc.velocity.v_per./particle.guidingcenter.gc.velocity.v;
particle.guidingcenter.gc.pitch_angle = asin(particle.guidingcenter.gc.lambda);
particle.guidingcenter.gc.larmor_radius = (particle.mass/abs(particle.charge))*particle.guidingcenter.gc.velocity.v_per./particle.guidingcenter.gc.B;
particle.guidingcenter.gc.angular_velocity = particle.charge*particle.guidingcenter.gc.B/particle.mass;
particle.guidingcenter.gc.angular_frequency = particle.guidingcenter.gc.angular_velocity/(2*pi);
particle.guidingcenter.gc.energy.kinetic_par = 0.5*particle.mass*particle.guidingcenter.gc.velocity.v_par.^2;
particle.guidingcenter.gc.energy.kinetic_per = 0.5*particle.mass*particle.guidingcenter.gc.velocity.v_per.^2;
particle.guidingcenter.gc.energy.kinetic = particle.guidingcenter.gc.energy.kinetic_par + particle.guidingcenter.gc.energy.kinetic_per;


% Evaluates the position of the guiding centre from the GC solution on the time array of 
% the full orbit solution in order to show the difference of the two
particle.error_x =  particle.orbit.guidingcenter_x - interp1(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.x, particle.orbit.t');
particle.error_y =  particle.orbit.guidingcenter_y - interp1(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.y, particle.orbit.t');
particle.error_z =  particle.orbit.guidingcenter_z - interp1(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.z, particle.orbit.t');


printf('%f s\n', toc());

