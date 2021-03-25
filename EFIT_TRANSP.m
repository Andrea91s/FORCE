%compare EFIT EQDSK with TRANSP EQDSK and the INTERNAL ONE FROM TRANSP
%{
TRANSP = load('TRANSP_EQDSK_29881_0245.mat');
EFIT = load('EFIT_EQDSK_29881_0245.mat');

addpath('/home/andrea/TRANSP')
filename  = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29881/29881O07.CDF';
TRANSP_BPOL =  nc_read(filename, 'BPOL');
TRANSP_PSI =  nc_read(filename, 'PLFLX');
TRANSP_R = nc_read(filename,'RZON');
TRANSP_Z = nc_read(filename,'LPOL');
%}
close all;

figure(1)
subplot(1,2,1)
pcolor(EFIT.EQD.R_grid, EFIT.EQD.Z_grid, EFIT.EQD.psi_2D)
hold all;
plot(EFIT.EQD.LCFS_R, EFIT.EQD.LCFS_Z, 'r', 'linewidth', 2);
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('Stream Function \Psi EFIT EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
%caxis([-0.15 0.06])


subplot(1,2,2)
pcolor(TRANSP.EQD.R, TRANSP.EQD.Z, TRANSP.EQD.psi_2D)
hold all;
plot(TRANSP.EQD.LCFS_R, TRANSP.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('Stream Function \Psi TRANSP EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
%caxis([-0.25 0.06])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(1,2,1)
pcolor(EFIT.EQD.R_grid, EFIT.EQD.Z_grid, EFIT.EQD.B.POL.R)
hold on;
plot(EFIT.EQD.LCFS_R, EFIT.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\theta,R} EFIT EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-0.16 0.16])

subplot(1,2,2)
pcolor(TRANSP.EQD.R, TRANSP.EQD.Z, TRANSP.EQD.B.POL.R)
hold on;
plot(TRANSP.EQD.LCFS_R, TRANSP.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\theta,R} TRANSP EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-0.16 0.16])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
subplot(1,2,1)
pcolor(EFIT.EQD.R_grid, EFIT.EQD.Z_grid, EFIT.EQD.B.POL.Z)
hold on;
plot(EFIT.EQD.LCFS_R, EFIT.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\theta,Z} EFIT EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-1 0.31])

subplot(1,2,2)
pcolor(TRANSP.EQD.R, TRANSP.EQD.Z, TRANSP.EQD.B.POL.Z)
hold on;
plot(TRANSP.EQD.LCFS_R, TRANSP.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\theta,Z} TRANSP EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-1 0.31])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
subplot(1,2,1)
pcolor(EFIT.EQD.R_grid, EFIT.EQD.Z_grid, EFIT.EQD.B.TOR)
hold on;
plot(EFIT.EQD.LCFS_R, EFIT.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\phi} EFIT EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-2 0])

subplot(1,2,2)
pcolor(TRANSP.EQD.R, TRANSP.EQD.Z, TRANSP.EQD.B.TOR)
hold on;
plot(TRANSP.EQD.LCFS_R, TRANSP.EQD.LCFS_Z, 'r', 'linewidth', 2); 
shading interp
colorbar
xlabel('R (m)')
ylabel('Z (m)')
title('B_{\phi} TRANSP EQDSK', 'fontsize', 12, 'fontweight', 'bold')
axis equal
xlim([0.2 1.75])
ylim([-1.54 1.54])
caxis([-2 0])

plot(EFIT.EQD.psi_rho_pol, EFIT.EQD.pressure, 'b', 'linewidth', 2,...
TRANSP.EQD.psi_rho_pol, TRANSP.EQD.pressure, 'r', 'linewidth', 2);
xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
ylabel('Pressure (Newton/m2)')
title('Plasma pressure', 'fontsize', 12, 'fontweight', 'bold')
legend('EFIT','TRANSP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
subplot(2,2,1)
plot(EFIT.EQD.R, EFIT.EQD.psi_1D, 'b', 'linewidth', 2,...
TRANSP.EQD.R, TRANSP.EQD.psi_1D, 'r', 'linewidth', 2);
xlabel('R (m)')
ylabel('stream function \Psi (Weber/rad)') 
title('Radial stream function \Psi', 'fontsize', 12, 'fontweight', 'bold')
legend('EFIT','TRANSP', 'location', 'southwest')



subplot(2,2,2)
plot(EFIT.EQD.psi_rho_pol, EFIT.EQD.F, 'b', 'linewidth', 2,...
TRANSP.EQD.psi_rho_pol, TRANSP.EQD.F, 'r', 'linewidth', 2);
xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
ylabel('Poloidal current function F (m Tesla)')
title('Poloidal current function', 'fontsize', 12, 'fontweight', 'bold')
legend('EFIT','TRANSP', 'location', 'southeast')


subplot(2,2,3)
plot(EFIT.EQD.psi_rho_pol, EFIT.EQD.q, 'b', 'linewidth', 2,...
TRANSP.EQD.psi_rho_pol, TRANSP.EQD.q, 'r', 'linewidth', 2);
xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
ylabel('Safety factor')
title('q safety factor', 'fontsize', 12, 'fontweight', 'bold')
legend('EFIT','TRANSP', 'location', 'southeast')


subplot(2,2,4)
plot(EFIT.EQD.R, EFIT.EQD.B.POL.R_radial, 'b--', 'linewidth', 2,...
     TRANSP.EQD.R, TRANSP.EQD.B.POL.R_radial, 'b', 'linewidth', 2,...
     EFIT.EQD.R, EFIT.EQD.B.POL.Z_radial, 'r--', 'linewidth', 2,...
     TRANSP.EQD.R, TRANSP.EQD.B.POL.Z_radial, 'r', 'linewidth', 2,...
     EFIT.EQD.R, EFIT.EQD.B.TOR_radial, 'g--', 'linewidth', 2,...
     TRANSP.EQD.R, TRANSP.EQD.B.TOR_radial, 'g', 'linewidth', 2);
xlabel('R (m)')
ylabel('Magnetic field (T)')
title('Radial profiles', 'fontsize', 12, 'fontweight', 'bold')
h=legend('EFIT B_{\theta,R} ','TRANSP B_{\theta,R}',...
       'EFIT B_{\theta,Z} ','TRANSP B_{\theta,Z}',...
       'EFIT B_{\phi} ','TRANSP B_{\phi}');
legend(h, 'location', 'eastoutside',"orientation", "vertical");
legend(h, 'boxoff');
set(h, 'fontsize', 9)     
axis([0.2 2 -1 1])

