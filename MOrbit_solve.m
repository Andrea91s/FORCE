function [particle] = MOrbit_solve(particle)
% solve the equation of motion

% ------------------------------------------------------------------------
% LSODE options
% ------------------------------------------------------------------------
set_lsode_options = 1;		% 0 for default (stiff)
                            % 1 for non-stiff
if (set_lsode_options == 1)
  lsode_options ('absolute tolerance', 1.0E-08); % default: 1.49012e-08
  lsode_options ('relative tolerance', 1.0E-08); % default: 1.49012e-08
  lsode_options('maximum order', -1);
  lsode_options('integration method', 'non-stiff');
  lsode_options('initial step size', -1);
  lsode_options('maximum step size', -1);
  lsode_options('minimum step size', 0);
  lsode_options('step limit', 100000);
elseif (set_lsode_options == 0)
  lsode_options ('absolute tolerance', 1.0e-08);
  lsode_options ('relative tolerance', 1.0e-08);
  lsode_options('maximum order', 2);
  lsode_options('integration method', 'stiff');
  lsode_options('initial step size', -1);
  lsode_options('maximum step size', -1);
  lsode_options('minimum step size', 0);
  lsode_options('step limit', 100000);
endif

% Store the info on the solver used
particle.orbit.solver.method = 'lsode';
particle.orbit.solver.abs_tol = lsode_options ('absolute tolerance');
particle.orbit.solver.rel_tol = lsode_options ('relative tolerance');
particle.orbit.solver.max_order = lsode_options('maximum order');
particle.orbit.solver.int_method =lsode_options('integration method');
particle.orbit.solver.init_step_size = lsode_options('initial step size');
particle.orbit.solver.max_step_size = lsode_options('maximum step size');
particle.orbit.solver.min_step_size = lsode_options('minimum step size');
particle.orbit.solver.step_limit =  lsode_options('step limit');

% ------------------------------------------------------------------------
% Defines global parameters
% ------------------------------------------------------------------------
global m = particle.mass;
global q = particle.charge; 

% ------------------------------------------------------------------------
% LSODE solution
% ------------------------------------------------------------------------
ic = [particle.initialconditions.x0; particle.initialconditions.y0; particle.initialconditions.z0; ...
      particle.initialconditions.vx0; particle.initialconditions.vy0; particle.initialconditions.vz0];

printf('%s', 'Full orbit solver: '); fflush(stdout());
tic();
[S, istate, msg] = lsode ("Lorentz_force", ic, particle.orbit.t);

% Compute solutions
particle.orbit.position.x = S(:,1);
particle.orbit.position.y = S(:,2);
particle.orbit.position.z = S(:,3);
particle.orbit.velocity.vx = S(:,4);
particle.orbit.velocity.vy = S(:,5);
particle.orbit.velocity.vz = S(:,6);


% ISTATE and Messages returned by LSODE
particle.orbit.ISTATE = istate;
particle.orbit.msg = msg;

particle.orbit.cpu_time = toc ();
printf('%f s\n', particle.orbit.cpu_time); fflush(stdout());
