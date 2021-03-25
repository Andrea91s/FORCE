function [rr, zz, ni3d, ti3d, r, ni, ti] = nc_D_density_temperature_plot(TIME, R, Z, NI3D, TI3D, time, fs, rho)
% Function that return the neutron emissivity for a given time and plots it
%
% INPUT
% TIME            time array [s]
% R               time x zones R x poloidal angle [cm]
% Z               time x zones Z x poloidal angle [cm]
% TI3D		array 		ion temperature time x length(R) x length(theta) coordinates in eV
% NI3D		array 		ion density time x length(R) x length(theta) coordinates in cm-3
% time            selected time for the plot [s]
%
%
% OUTPUT
% rr			      R at selected time [cm]
% zz			      Z at selected time [cm]
% ni3d			    ion density length(R) x length(theta) coordinates in cm-3
% ti3d          ion temperature length(R) x length(theta) coordinates in eV
% r             radial profile major radius at Z = 0 [cm]
%
%
% Example:
% filename = '29880U32.CDF';
% directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U32/';
% [TIME, TI3D, NI3D, R, Z, RAXIS, ZAXIS, UTI3D, UNI3D] = nc_D_density_temperature(strcat(directory,filename));
% [rr, zz, ni3d, ti3d, r, ni, ti] = nc_D_density_temperature_plot(TIME, R, Z, NI3D, TI3D, time = 0.26, fs = 1, rho);

% Find the index corresponding to the selected time
%idx = max(find(TIME <= time))
idx=56

% Select the time
rr = squeeze(R(idx,:,:));
zz = squeeze(Z(idx,:,:));
ni3d = squeeze(NI3D(idx,:,:));
ti3d = squeeze(TI3D(idx,:,:));


% Calculate the radial profile for Z = 0 
r = linspace(min(rr(:)), max(rr(:)), 1000);
ti = griddata(rr,zz,ti3d,r,0);
ni = griddata(rr,zz,ni3d,r,0);


% Plot the results
%load MASTcmap2.dat

% Density 2D map
figure(10)
pcolor(rr, zz, ni3d)
shading flat
if (fs == 1)
    hold on
	for k = 1:5:60
        plot(rr(k,:), zz(k,:), 'k')
    end
	hold off
endif

%set(gcf,'Colormap',MASTcmap2) 
xlabel('R (cm)', 'fontsize', 12)
ylabel('Z (cm)', 'fontsize', 12)
axis equal
axis([0 200 -150 150])
h = colorbar;
set(h, 'title', 'n_D [cm^{-3}]')
title (['D density at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)

% Temperature 2D map
figure(11)
pcolor(rr, zz, ti3d)
shading flat
if (fs == 1)
    hold on
	for k = 1:5:60
        plot(rr(k,:), zz(k,:), 'k')
    end
	hold off
endif
rho_n =squeeze(ni3d(:,1));
rho_t =squeeze(ti3d(:,1));
%set(gcf,'Colormap',MASTcmap2) 
xlabel('R (cm)', 'fontsize', 12)
ylabel('Z (cm)', 'fontsize', 12)
axis equal
axis([0 200 -150 150])
h = colorbar;
set(h, 'title', 'n_D [eV]')
title (['D temperature at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)

figure(12)
plot(r, ti, 'k-', 'linewidth', 2);
xlabel('R (cm)', 'fontsize', 12)
ylabel('T_D (eV)', 'fontsize', 12)
title (['D temperature at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)


figure(13)
plot(r, ni, 'k-', 'linewidth', 2);
xlabel('R (cm)', 'fontsize', 12)
ylabel('n_D [cm^{-3}]', 'fontsize', 12)
title (['D density at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)

figure(14)
plot(rho, rho_n, 'k-', 'linewidth', 2);
xlabel('\rho ', 'fontsize', 12)
ylabel('n_D [cm^{-3}]', 'fontsize', 12)
title (['D density at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)


figure(15)
plot(rho, rho_t, 'k-', 'linewidth', 2);
xlabel('\rho ', 'fontsize', 12)
ylabel('T [eV]', 'fontsize', 12)
title (['temperature at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)