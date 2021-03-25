function [particle] = MOrbit(PTCL, POS, VEL_or_EPA, TIME, EQLST, SGC, plotyn, saveyn)
saveyn=1;
% Orbit following code
% All unit in MKS apart from energy (eV)
% 
% INPUT
% PTCL		    particle properties: -1 electron, 1 = proton, 2 = deuteron, 3 = tritium
% PSO		    particle initial position vector [x0 y0 z0]
% VEL_or_EPA	particle initial:  
%               1] ROW VECTOR = velocity vector [vx0 vy0 vz0] in m/s 
%               2] COL VECTOR = energy , pitch angle, gyro-angle [E; alpha; phi] in (eV, deg, rad)
%               NB For a random gyro-angle set phi = 2*pi*rand(1);
% EQLST   	    2nd output of the efit_read.m routine with EFIT data for the given time slice
% TIME		    LSode time vector [t_start t_end step]
% SGC		    keyword = 0 does NOT solve using guiding centre approximation
%                         1 solve guiding cente equation of motion (approximation)
% plotyn        % 0 no plots generated
%                 1 plots generated
%                 2 plots generated on the external monitor assumed to be on the right hand side of a
%                   full HD monitor
%
% EXAMPLE:
%
% 1) Read EFIT equilibrium from HDF5 file directly or from a pre-read and stored file
%    [EQL, EQLST] = efit_read(filename = '~/Documents/MAST/EFIT/29881/efitOut.hdf5', st = 0.25, plotyn = 0, saveyn = 0);
%    clear -x EQLST
%    save EFIT_29881_at_0.250_s.mat
%    load EFIT_29881_at_0.250_s.mat
%
% 2) Alternatively used IDAM EFIT equilibria
%    [EQL, EQLST] = EFIT_IDAM(shot = 29881, ts = 0.245, off_on_line = 0, plotyn = 0, n_lvl = 0, min_lev = 0);
%    clear -x EQLST
% 
% 3) Reading from EQDSK files
%    filename = '/home/marco/Documents/MAST/EFIT/29976/U39/g29976.00210';
%    EQD = read_eqdsk(filename, 0);
%    EQLST = EQD2EQLST(EQD);
%    clear EQD filename
%
% 4) Run the code:
%    Case 1: Deuteron 1 keV, pitch angle = 0 DEG (parallel to B); random gyrophase
%    [particle] = MOrbit(PTCL = 2, POS = [1.2 0 0], VEL_or_EPA = [1000; 0; 2*pi*rand(1)], TIME = [0 1E-5 5000], EQLST, SGC = 0 , plotyn = 1, saveyn = 0);
%
%    Case 2: Deuteron 60 keV, pitch angle = 180 DEG (anti-parallel to B); random gyrophase
%    [particle] = MOrbit(PTCL = 2, POS = [1.0 0 0], VEL_or_EPA = [60000; 180; 0], TIME = [0 1E-5 5000], EQLST, SGC = 0 , plotyn = 1, saveyn = 0);
%
%    Case 3: Deuteron 60 keV, pitch angle = 0 DEG (parallel to B); random gyrophase
%    [particle] = MOrbit(PTCL = 2, POS = [1.0 0 0], VEL_or_EPA = [60000; 0; 0], TIME = [0 1E-5 5000], EQLST, SGC = 0 , plotyn = 1, saveyn = 0);


% Add the required path
addpath('~/Documents/Fusion/', '~/Documents/MAST/EFIT', '~/Documents/MAST/Octave', '~/Documents/MAST/LINE2/');



% ------------------------------------------------------------------------
% Set the EFIT magnetic variables as global and Physical constants
% ------------------------------------------------------------------------
MOrbit_global_variables;


% ------------------------------------------------------------------------
% Particle initial conditions
% ------------------------------------------------------------------------
[particle] = MOrbit_pre_processor(PTCL, POS, VEL_or_EPA, TIME, SGC);


% ---------------------------------------------------
% SOLVER
% ---------------------------------------------------
[particle] = MOrbit_solve(particle);


% ------------------------------------------------------------------------
% FULL ORBIT POST PROCESSOR
% ------------------------------------------------------------------------
[particle] = MOrbit_post_processor(particle);


% ------------------------------------------------------------------------
% GUIDING CENTRE CALCULATIONS from ORBIT calculations
% ------------------------------------------------------------------------
[particle] = MOrbit_guiding_centre(particle);


% ------------------------------------------------------------------------
% GYRO-AVERARGE calculations
% ------------------------------------------------------------------------
%[particle] = MOrbit_gyro_average(particle);


% ------------------------------------------------------------------------
% GUIDING CENTRE CALCULATIONS
% ------------------------------------------------------------------------
if (particle.guidingcenter.gc.solve == 1)
    [particle] = MOrbit_guiding_center_approx(particle);
endif


     
% ------------------------------------------------------------------------
% Turning points location, bounce frequency and mirror ratio
% ------------------------------------------------------------------------
%[particle] = MOrbit_turning_points(particle);


% ------------------------------------------------------------------------
% Toroidal velocity reflection or not?
% ------------------------------------------------------------------------
%[particle] = MOrbit_vtor_check(particle);


% ------------------------------------------------------------------------
% Confined or lost ?
% ------------------------------------------------------------------------
%[particle] = MOrbit_lost_particles(particle);


% ------------------------------------------------------------------------
% Encircling the magnetic axis or not?
% ------------------------------------------------------------------------
%[particle] = MOrbit_cyrcling_particles(particle, EQLST);


% ------------------------------------------------------------------------
% DRIFT/Precession frequency
% ------------------------------------------------------------------------
% This depends really on the problem at hands so maybe not worth trying a
% catch it all solution.

% ------------------------------------------------------------------------
% Orbit classification
% ------------------------------------------------------------------------
printout = 0;
[particle] = MOrbit_classification(particle, EQLST, printout);


% ------------------------------------------------------------------------
% PRINT THE OUTPUTS
% ------------------------------------------------------------------------
if (plotyn > 0)
   % MOrbit_printout(particle)
   MOrbit_printout(particle, EQLST, sfb = 0);
endif
    


% ------------------------------------------------------------------------
% PLOT THE OUTPUTS
% ------------------------------------------------------------------------
% Plot the orbit
if (plotyn > 0)
  MOrbit_plot(particle, EQLST, plotyn);
endif
    


% ------------------------------------------------------------------------
% SAVE THE OUTPUT
% ------------------------------------------------------------------------
if (saveyn == 1)
  MOrbit_printout(particle, EQLST, sfb = 1);
  MOrbit_save(particle, EQLST);
endif




































