function [FieldLine] = Bfield_Line_Solver(P, L, EQLST)
% solve the magnetic field line equation
%
% Inputs
% P         initial position
% L         field line length array;
%
%
% Example;
% L = linspace (0, 100, 5000);                          % 100 m long field line, 5000 steps
% P = [1.4; 0; 0];                                      % starting from this point in [m]
% load EFIT_29881_at_0.250_s.mat;                       % load the magnetic field equilibirum
% MOrbit_global_variables;                              % declare the global MOrbit_global_variables
% [FieldLine] = Bfield_Line_Solver(P, L, EQLST);        % solve & plot the B field line equations

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
  lsode_options ('absolute tolerance', 1.49012e-08);
  lsode_options ('relative tolerance', 1.49012e-08);
  lsode_options('maximum order', 2);
  lsode_options('integration method', 'stiff');
  lsode_options('initial step size', -1);
  lsode_options('maximum step size', -1);
  lsode_options('minimum step size', 0);
  lsode_options('step limit', 100000);
endif

% Store the info on the solver used
FieldLine.solver.method = 'lsode';
FieldLine.solver.abs_tol = lsode_options ('absolute tolerance');
FieldLine.solver.rel_tol = lsode_options ('relative tolerance');
FieldLine.solver.max_order = lsode_options('maximum order');
FieldLine.solver.int_method =lsode_options('integration method');
FieldLine.solver.init_step_size = lsode_options('initial step size');
FieldLine.solver.max_step_size = lsode_options('maximum step size');
FieldLine.solver.min_step_size = lsode_options('minimum step size');
FieldLine.solver.step_limit =  lsode_options('step limit');

% ------------------------------------------------------------------------
% Defines global parameters
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% LSODE solution
% ------------------------------------------------------------------------
FieldLine.IC = P;
FieldLine.L = L;

printf('%s', 'Magnetic field line solver: '); fflush(stdout());
tic();
[S, istate, msg] = lsode ("Bfield_Line", FieldLine.IC, FieldLine.L);

% Compute solutions
FieldLine.x = S(:,1);
FieldLine.y = S(:,2);
FieldLine.z = S(:,3);


% ISTATE and Messages returned by LSODE
FieldLine.ISTATE = istate;
FieldLine.msg = msg;

FieldLine.cpu_time = toc ();
printf('%f s\n', FieldLine.cpu_time); fflush(stdout());
printf('%s \n', FieldLine.msg); fflush(stdout());


% Plot
MAST_Geometry_3D(EQLST);
hold on
plot3(S(:,1), S(:,2), S(:,3),'r', 'linewidth', 2)
hold off

