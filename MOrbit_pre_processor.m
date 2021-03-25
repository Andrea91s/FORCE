function [particle] = MOrbit_pre_processor(PTCL, POS, VEL_or_EPA, TIME, SGC)
% ------------------------------------------------------------------------
% Particle initial conditions
% ------------------------------------------------------------------------

tic();
printf('%s', 'Full orbit pre-processor: '); fflush(stdout());



% ------------------------------------------------------------------------
% Set keyword for the guiding center approximation
% ------------------------------------------------------------------------
particle.guidingcenter.gc.solve = SGC;	% 0 = No
                                        % 1 = Yes
                                        % If YES then the time step for solving the
                                        %        GC equation is taken as 1/10 of that
                                        %        of the FO solver.

% ------------------------------------------------------------------------
% Particle initial conditions
% ------------------------------------------------------------------------


% Initial position in carthesian coordiantes (unit in m)
particle.initialconditions.x0 = POS(1);
particle.initialconditions.y0 = POS(2);
particle.initialconditions.z0 = POS(3);


% Solve the equation of motion for time (in sec)
particle.times.start_time = TIME(1);
particle.times.end_time = TIME(2); %0.5*4e-5;
particle.times.steps = TIME(3);   %5000;

% Simulation time step:
% 1. entered directly as:
%    particle.times.dt = 1E-9;
%
% 2. depending on the initial gyration frequency and on the number of time points
%    per orbit
%     particle.times.steps_per_gyration = 300;
%     particle.times.dt = (1/particle.initialconditions.frequency)/particle.times.steps_per_gyration;

% Number of time steps in simulation
%particle.times.steps = round((particle.times.end_time - particle.times.start_time)/particle.times.dt);

% Simulation time
particle.orbit.t = linspace(particle.times.start_time, particle.times.end_time, particle.times.steps)';
particle.times.dt = (particle.times.end_time - particle.times.start_time)/particle.times.steps;


% The keyword below specify is (vx0, vy0, vz0) or (E0, alpha0)
% should be used to define the particle initial conditions:
[rown, coln] = size(VEL_or_EPA);
switch rown
  case 3
      particle.initialconditions.choose_initial_velocity = 1;
      particle.initialconditions.energy_eV = VEL_or_EPA(1);
      particle.initialconditions.pitch_angle = VEL_or_EPA(2)*(pi/180);
      particle.initialconditions.phi = VEL_or_EPA(3); 
  case 1
      particle.initialconditions.choose_initial_velocity = 0;
      particle.initialconditions.vx0 = VEL_or_EPA(1);
      particle.initialconditions.vy0 = VEL_or_EPA(2);
      particle.initialconditions.vz0 = VEL_or_EPA(3);
endswitch   


% ------------------------------------------------------------------------
% % Particle properties
% ------------------------------------------------------------------------
switch PTCL
  case -1
    particle.mass = physical_constant ('electron mass');
    particle.charge = -physical_constant ('elementary charge');
  case 1
    particle.mass = physical_constant ('proton mass');
    particle.charge = physical_constant ('elementary charge');
  case 2
    particle.mass = physical_constant ('deuteron mass');
    particle.charge = physical_constant ('elementary charge');  
  case 3
    particle.mass = 3*physical_constant ('triton mass');
    particle.charge = physical_constant ('elementary charge');
endswitch


% Calculate the magnetic field versors at the initial position
particle.initialconditions.Bx = Bfield(particle.initialconditions.x0, particle.initialconditions.y0, particle.initialconditions.z0)(1);
particle.initialconditions.By = Bfield(particle.initialconditions.x0, particle.initialconditions.y0, particle.initialconditions.z0)(2);
particle.initialconditions.Bz = Bfield(particle.initialconditions.x0, particle.initialconditions.y0, particle.initialconditions.z0)(3);
particle.initialconditions.B = sqrt(particle.initialconditions.Bx^2 + particle.initialconditions.By^2 + particle.initialconditions.Bz^2);
particle.initialconditions.bx =  particle.initialconditions.Bx /particle.initialconditions.B; 
particle.initialconditions.by =  particle.initialconditions.By /particle.initialconditions.B; 
particle.initialconditions.bz =  particle.initialconditions.Bz /particle.initialconditions.B; 
particle.initialconditions.d = sqrt(particle.initialconditions.bx^2 + particle.initialconditions.by^2);

% Calculate the electric field versors at the initial position
[particle.initialconditions.Ex, particle.initialconditions.Ey, particle.initialconditions.Ez, particle.initialconditions.E] = Efield(particle.initialconditions.x0, particle.initialconditions.y0, particle.initialconditions.z0);
if (particle.initialconditions.E != 0)
    particle.initialconditions.Ex =  particle.initialconditions.Ex/particle.initialconditions.E; 
    particle.initialconditions.Ey =  particle.initialconditions.Ey/particle.initialconditions.E; 
    particle.initialconditions.Ez =  particle.initialconditions.Ez/particle.initialconditions.E; 
else
    particle.initialconditions.Ex =  0; 
    particle.initialconditions.Ey =  0; 
    particle.initialconditions.Ez =  0; 
end


% Particle initial gyration frequency and angular frequency:
particle.initialconditions.omega = particle.charge*particle.initialconditions.B/particle.mass;
particle.initialconditions.frequency = particle.initialconditions.omega/(2*pi);

% Initial velocity: this can be set by the user or calculated from energy and pitch-angle values
% by setting the following keyowrd to 0 or 1 respectively

if (particle.initialconditions.choose_initial_velocity == 0)
    % user input initial velocity
    particle.initialconditions.v = sqrt(particle.initialconditions.vx0^2 + particle.initialconditions.vy0^2+particle.initialconditions.vz0^2);
    particle.initialconditions.energy_J = 0.5*particle.mass*particle.initialconditions.v^2;
    particle.initialconditions.energy_eV =  energy_unit (particle.initialconditions.energy_J, conversion = 0);
    [ParComp, PerComp, Angle] = projection([particle.initialconditions.vx0 particle.initialconditions.vy0 particle.initialconditions.vz0], [particle.initialconditions.Bx particle.initialconditions.By particle.initialconditions.Bz]);
    particle.initialconditions.v_pll_x = ParComp(1);
    particle.initialconditions.v_pll_y = ParComp(2);
    particle.initialconditions.v_pll_z = ParComp(3);
    particle.initialconditions.v_per_x = PerComp(1);
    particle.initialconditions.v_per_y = PerComp(2);
    particle.initialconditions.v_per_z = PerComp(3);    
    particle.initialconditions.v_parl = norm(ParComp);
    particle.initialconditions.v_perp = norm(PerComp);  
    particle.initialconditions.pitch_angle = Angle;
    particle.initialconditions.lambda = cos(particle.initialconditions.pitch_angle);
    
    % NB! The gyration angle below is the one obtained in the ref. frame in which:
    %     X' = v x B
    %     Y' = B x (v x B)
    %     Z' = B
    particle.initialconditions.phi = pi/2;

elseif (particle.initialconditions.choose_initial_velocity == 1)
    % Energy and pitch angle initial values
    particle.initialconditions.energy_J = energy_unit (particle.initialconditions.energy_eV, conversion = 1);
    particle.initialconditions.v = sqrt(2*particle.initialconditions.energy_J/particle.mass);
    particle.initialconditions.v_perp =  particle.initialconditions.v*sin(particle.initialconditions.pitch_angle);
    particle.initialconditions.v_parl =  particle.initialconditions.v*cos(particle.initialconditions.pitch_angle);
    particle.initialconditions.lambda = cos(particle.initialconditions.pitch_angle);
   
    % Calculate the inital velocity components from Energy and pitch angle
    if (particle.initialconditions.bx == 0 && particle.initialconditions.by == 0 && particle.initialconditions.bz != 0)
      particle.initialconditions.vx0 = particle.initialconditions.v_perp*cos(particle.initialconditions.phi);
      particle.initialconditions.vy0 = particle.initialconditions.v_perp*sin(particle.initialconditions.phi);
      particle.initialconditions.vz0 = particle.initialconditions.v_parl;
    else
      particle.initialconditions.vx0 = (particle.initialconditions.by/particle.initialconditions.d)*particle.initialconditions.v_perp*cos(particle.initialconditions.phi) + ...
				      (particle.initialconditions.bx*particle.initialconditions.bz/particle.initialconditions.d)*particle.initialconditions.v_perp*sin(particle.initialconditions.phi) + ...
				      particle.initialconditions.bx*particle.initialconditions.v_parl;
      particle.initialconditions.vy0 = -(particle.initialconditions.bx/particle.initialconditions.d)*particle.initialconditions.v_perp*cos(particle.initialconditions.phi) + ...
				      (particle.initialconditions.by*particle.initialconditions.bz/particle.initialconditions.d)*particle.initialconditions.v_perp*sin(particle.initialconditions.phi) + ...
				      particle.initialconditions.by*particle.initialconditions.v_parl;
      particle.initialconditions.vz0 = - particle.initialconditions.d*particle.initialconditions.v_perp*sin(particle.initialconditions.phi) +  particle.initialconditions.bz*particle.initialconditions.v_parl;                
    endif

endif

% Initial larmor radius
particle.initialconditions.larmor_radius = (particle.mass/abs(particle.charge))*particle.initialconditions.v_perp./particle.initialconditions.B;

% NB! The gyration angle below is the one obtained in the ref. frame in which:
%     X' = v x B
%     Y' = B x (v x B)
%     Z' = B
% particle.initialconditions.phi = atan2(particle.initialconditions.vy0, particle.initialconditions.vx0);

printf('%f s\n', toc()); fflush(stdout());

