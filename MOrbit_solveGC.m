function [particle] = MOrbit_solveGC(particle)

% ------------------------------------------------------------------------
% LSODE options
% ------------------------------------------------------------------------
set_lsode_options = 1;		% 0 for default (stiff)
                            % 1 for non-stiff
if (set_lsode_options == 1)
  lsode_options ('absolute tolerance', 1.0E-06); % default: 1.49012e-08
  lsode_options ('relative tolerance', 1.0E-06); % default: 1.49012e-08
  lsode_options('maximum order', -1);
  lsode_options('integration method', 'non-stiff');
  lsode_options('initial step size', -1);
  lsode_options('maximum step size', -1);
  lsode_options('minimum step size', 0);
  lsode_options('step limit', 100000);
elseif (set_lsode_options == 0)
  lsode_options ('absolute tolerance', 1.49012e-08);
  lsode_options ('relative tolerance', 1.49012e-08);
  lsode_options('maximum order', 5);
  lsode_options('integration method', 'stiff');
  lsode_options('initial step size', -1);
  lsode_options('maximum step size', -1);
  lsode_options('minimum step size', 0);
  lsode_options('step limit', 100000);
endif

% Store the info on the solver used
particle.guidingcenter.solver.method = 'lsode';
particle.guidingcenter.solver.abs_tol = lsode_options ('absolute tolerance');
particle.guidingcenter.solver.rel_tol = lsode_options ('relative tolerance');
particle.guidingcenter.solver.max_order = lsode_options('maximum order');
particle.guidingcenter.solver.int_method =lsode_options('integration method');
particle.guidingcenter.solver.init_step_size = lsode_options('initial step size');
particle.guidingcenter.solver.max_step_size = lsode_options('maximum step size');
particle.guidingcenter.solver.min_step_size = lsode_options('minimum step size');
particle.guidingcenter.solver.step_limit =  lsode_options('step limit');

% Particle velocity (without E field this is conserved)


ic = [particle.guidingcenter.gc.position.x0; particle.guidingcenter.gc.position.y0; particle.guidingcenter.gc.position.z0; ...
      particle.guidingcenter.gc.position.v_par0];  
      
particle.guidingcenter.gc.t = linspace(particle.guidingcenter.gc.start_time, particle.guidingcenter.gc.end_time, particle.guidingcenter.gc.steps);


% ------------------------------------------------------------------------
% LSODE solution
% ------------------------------------------------------------------------    
particle.guidingcenter.gc.time_steps = int32(t_step/10);
particle.guidingcenter.gc.t = linspace(particle.times.start_time, particle.times.end_time, particle.guidingcenter.gc.time_steps);
tic();
printf('%s', 'Guiding cente solver: ');
S = lsode ("Lorentz_forceGC", ic, particle.guidingcenter.gc.t);


% Compute solutions

particle.guidingcenter.gc.position.x = S(:,1);
particle.guidingcenter.gc.position.y = S(:,2);
particle.guidingcenter.gc.position.z = S(:,3);
particle.guidingcenter.gc.velocity.v_par = S(:,4)';
clear S

particle.guidingcenter.cpu_time = toc ();
printf('%f s\n', particle.guidingcenter.cpu_time);
