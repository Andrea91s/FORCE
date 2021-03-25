function [rho, TIME, TI3D, NI3D, R, Z, RAXIS, ZAXIS, UTI3D, UNI3D, PFLUX] = nc_D_density_temperature(filename)
% Script to access Deuterium density and temperature from TRANSP output file
% INPUT
% filename	string		name of the NETCDF file containing TRANSP output
% plotyn 	integer		keyword: 0 no plot, 1 plot
%
% OUTPUT
% TIME		array		  time in sec
% TI3D		array 		ion temperature time x length(R) x length(theta) coordinates in eV
% NI3D		array 		ion density time x length(R) x length(theta) coordinates in cm-3
% R		    array 		flux surface radial coordinate in cm
% Z		    array 		flux surface vertical coordinate in cm
% RAXIS		array 		magnetic axis radial coordinate in cm
% ZAXIS 	array 		magnetic axis vertical coordinate in cm
% 
% Example:
% filename = '29881O07.CDF';
% filename = '99999I20.CDF';
% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/29909C18.CDF';
% [rho, TIME, TI3D, NI3D, R, Z, RAXIS, ZAXIS, UTI3D, UNI3D, PFLUX] = nc_D_density_temperature(filename);


% Read the time
TIME = nc_read(filename, 'TIME');   % time in sec
X = ncread(filename, 'X');
rho = squeeze(X(:,1));
% Specify how to create the 2D neutron emissivity
NTHETA = 50;
method = 1;

% Read the total neutron emissivities defined on the centres of
% the flux surfaces boundaries
[TI3D, UTI3D] = nc_neutronemissivity(filename, 'TI', NTHETA, method);
[NI3D, UNI3D] = nc_neutronemissivity(filename, 'ND', NTHETA, method);

% Read the flux surfaces boundaries, the normalized boundaries and centres: 
%   NOTE: 1) R,Z are the coordinates of the orignal flux surfaces WITHOUT the
%         additionl flux surface at the magnetic axis which is needed to avoid
%         a "hole" in the 2D pcolor plot of the neutron emissivity
%
%         2) UR, UZ are the coordinates of the orignal flux surfaces WITH the
%         additionl flux surface at the magnetic axis which is needed to avoid
%         a "hole" in the 2D pcolor plot of the neutron emissivity

[R, Z, RAXIS, ZAXIS, PFLUX, UR, UZ] = nc_fluxsurfaces(filename);


% Note that the commands below should be matched to the method used to calculate
% the neutron emissivities
if (method == 0) 
    clear R Z
    R = UR;
    Z = UZ;
    clear UR UZ
    [m,n] = size(DV);
    DV = [zeros(m,1) DV];
endif















	
