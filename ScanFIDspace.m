% Scan the fast ion distribution function Energy - Pitch Angle space
% launching the particle from a given point in space

pos=[0.7];

% Select the particle
PTCL = 2;

for abc=1:length(pos)
% Pick the particle launching position POS = [X Y Z]
% with units of meters
POS = [pos(abc) 0 0];

% Particle Energies in keV
Energies = 1000*[50];
NE = length(Energies);

% Particle Pitch angle in deg
%lambda = [0.5];
%PA = acos(lambda)*180/pi;
PA = [60];
NP = length(PA);

% Particle gyro-angle
phi = [0];

% Integration time TIME = [tstart tend steps] in second
TIME = [0 1E-5 20000];

% No plots
SGC = 0;
plotyn = 0;
plotyn_text=0;

% But save the data please
saveyn = 1;


% Load EFIT data
load test.mat;

saveyn

% Loop over all energies and pitch angles

  for ne = 1:NE
    for np = 1:NP
      for gyro =1:1
        clear particle
        printf('Solving for E = %2.0f keV, PA = %d ...\n', Energies(ne)/1000, PA(np)); fflush(stdout());
        [particle] = MOrbit(PTCL, POS, VEL_or_EPA = [Energies(ne); PA(np); phi(gyro)], TIME, EQLST, SGC, plotyn, plotyn_text, saveyn);
        printf('... done!\n\n', Energies(ne)/1000, PA(np), phi(gyro)); fflush(stdout());
      endfor
    endfor
  endfor
endfor