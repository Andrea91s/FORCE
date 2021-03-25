% Script for the evaluation of the fractions of the different populations as
% a function of the gyro-angle csi

close all;
clear all;

% ----------------------------------------------------------------------------------------------------
% Select case
% ----------------------------------------------------------------------------------------------------
% Load a fast ion distribution function and calculates the fast ions
% populations using the old method (interval based)

% Name of the file from which the data should be read
transp_runs = [16, 31, 32];
idx = 3;

% Choose either pre or post
sawtooth_phase = 'pre';
%sawtooth_phase = 'post';

% Choose radial position (Z = 0 for all analyses here)
Radial_Position = [90 100 110 120];
idx_R = 3;


% Directory where the FID are stored
directory = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29880/';
switch (sawtooth_phase)
    case {'pre'}
        fid = 1;
    case {'post'}
        fid = 2;
endswitch
filename =  strcat('U', num2str(transp_runs(idx)), '/29880U', num2str(transp_runs(idx)), '_fi_', num2str(fid) , '.cdf');
filename =  strcat(directory, filename);

addpath('/home/andrea/TRANSP/');
% Read the data for the selected point in space

[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id] = fastions(filename, RR = Radial_Position(idx_R), RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
% Save the selected FI distribution for plotting in Veusz
w = squeeze(FID(id,:,:))';
w=flip(w,2);

% ----------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------
% Load the orbit classification matrix 
% ----------------------------------------------------------------------------------------------------
%filename = '/home/marco/Documents/Orbit/MASTOrbit/Results/29880_Orbit_Classification_R_090.mat';
%classification_matrix = load ('-ascii', filename, 'classification_matrix');

load /home/andrea/Documents/MASTOrbit/29880_Orbit_Classification_R_110.mat

% Defines E and Lambda used in the classification matrix
energy = linspace(60, 5, 12)*1000;      % energy [eV]
de = (5)*1000;

lambda = linspace(-1,0,11);             % pitch
dl = 0.1;

% Define csi
csi = 30:30:330;
NC = numel(csi);



% ----------------------------------------------------------------------------------------------------
% Iterate over CSI to calculate the fractions
% ----------------------------------------------------------------------------------------------------
BI = zeros(NC,5);

orbit_type = 1:5;
for nc = 1:NC

    % Calculate the densities for all types of orbits
    for k = 1:5
        [BI(nc,k)] = bin_integral(squeeze(classification_matrix(nc,:,:)), orbit_type(k), energy, de, lambda, dl, EN, PA, fid = w, check = 1);
    end
    
    % Calculate the total density of fast ions that have ab orbit topology of the chosen type
    BIT(nc) = sum(BI(nc,:));
    
    % Calculate the ratios
    for k = 1:5
        Ratios(nc,k) =  BI(nc,k)/BIT(nc);
    end
end

PassingRatios = sum(Ratios(:,1:2),2);
TrappedRatios = sum([Ratios(:,3) Ratios(:,5)],2);
% ----------------------------------------------------------------------------------------------------

% ----------------------------------------------------------------------------------------------------
% Plot the results
% ----------------------------------------------------------------------------------------------------
figure(100)
plot(csi, Ratios, {'or', 'og', 'ok', 'or', 'ob'}, 'linewidth', 2)
xlabel('gyro-phase angle (deg)', 'fontsize', 12)
legend('Passing', 'Stagnation', 'Banana', 'Lost', 'Potato')
ylabel('Fraction', 'fontsize', 12)
title(['Pulse 29880 at R =' num2str(Radial_Position(idx_R)) ' m'])

figure(200)
plot(csi, PassingRatios, 'ob', 'linewidth', 2, csi, TrappedRatios, 'og', 'linewidth', 2, csi, Ratios(:,4), 'or', 'linewidth', 2)
xlabel('gyro-phase angle (deg)', 'fontsize', 12)
legend('Passing', 'Trapped', 'Lost')
ylabel('Fraction', 'fontsize', 12)
title(['Pulse 29880 at R =' num2str(Radial_Position(idx_R)) ' m'])
csi=csi';
save('csi.dat','csi','-ascii')
save('frac110.dat','PassingRatios','-ascii')
return




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




