% Script testing the bin-based evaluation of the fast ion densities in different
% regions of the energy-pitch space


% Load a fast ion distribution function and calculates the fast ions
% populations using the old method (interval based)
fid_results;


% Load the orbit classification matrix 
directory = '/home/marco/Documents/Orbit/MASTOrbit/';
filename = 'Orbit_Classification_29880_at_0.255_s_R0.90_Z0.00_manual.mat';
load ('-ascii', strcat(directory, filename));
classification_matrix = Orbit_Classification_29880_at_0_255_s_R0_90_Z0_00_manual;
clear Orbit_Classification_29880_at_0_255_s_R0_90_Z0_00_manual

% Defines E and Lambda used in the classification matrix
energy = linspace(60, 5, 12)*1000;      % energy [eV]
de = (5)*1000;

lambda = linspace(-1,1,21);             % pitch
dl = 0.1;


[BIT, BI] = bin_integral(classification_matrix, orbit_type, energy, de, lambda, dl, EN, PA, fid, check) 

% Example

% Calculate the densities for all types of orbits
orbit_type = 1:5;
for k = 1:5
    [BI(k)] = bin_integral(classification_matrix, orbit_type(k), energy, de, lambda, dl, EN, PA, fid = w, check = 0);
end 

% Calculate the total density
[BIT] = bin_integral(classification_matrix, -1, energy, de, lambda, dl, EN, PA, fid = w, check = 0);

% Print results
fprintf('------------------------------------------------------\n')
fprintf('Densities from bin integrals\n')
fprintf('------------------------------------------------------\n')
fprintf('Total: %1.4g\n', BIT)
for k = 1:5
    fprintf('Orbit type %d: %1.4g (%1.2f)\n', orbit_type(k), BI(k), 100*BI(k)/BIT);
end 

% Strange difference on the absolute levels of a factor of approximate 10... could it be related
% to the difference in bin widths??
DE = diff(EN(1:2));     % TRANSP energy bin width
DL = diff(PA(1:2));     % TRANSP pitch bin width

bin_size_ratio = de*dl/(DE*DL);

% Print results corrected by the bin_size ratio
fprintf('------------------------------------------------------\n')
fprintf('Densities from bin integrals\n')
fprintf('------------------------------------------------------\n')
fprintf('Bin size ratio: %.14f\n', bin_size_ratio)
fprintf('Total: %1.4g\n', BIT/bin_size_ratio)
for k = 1:5
    fprintf('Orbit type %d: %1.4g (%1.2f)\n', orbit_type(k), BI(k)/bin_size_ratio, 100*BI(k)/BIT);
end 




