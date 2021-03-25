function [] = plot_eqdsk(EQD, filename);
% function that plots the equilibrium parameters from 
% an EQDSK file
%
% INPUT
% EQD 		structure containing the EFIT equilibria parameters
% filename	name of the file containing the equilibirum data




% ----------------------------------------------------------------------------------------------------
% PLOT G-FILES quantities
% ----------------------------------------------------------------------------------------------------

figure (1, 'position', [100 100 1400 800])
  
  % Plot the stream function PSI
  subplot(2,3,1)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the stream function
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.psi_2D); 
    title('Stream Function \Psi', 'fontsize', 12, 'fontweight', 'bold')
    shading interp
    h = colorbar;
    set(h, 'title', '(Weber/rad)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, 30, 'k', 'linewidth', 1.);
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal

    hold off

  
  % Plot the normalized stream function PSI
  subplot(2,3,2)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the stream function
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.psi_2D_norm); 
    title('Normalized Stream Function \Psi_{N}', 'fontsize', 12, 'fontweight', 'bold')
    shading interp
    h = colorbar;

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, 30, 'k', 'linewidth', 1.);
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    
    
    
  % Plot the poloidal flux PSI [Weber]
  subplot(2,3,3)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the poloidal flux
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D); 
    shading interp
    title('Poloidal Flux \Psi_{p}', 'fontsize', 12, 'fontweight', 'bold')
    h = colorbar;
    set(h, 'title', '(Weber)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, 30, 'k', 'linewidth', 1.);
    contour(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal

    hold off
      
    
  % Plasma pressure
  subplot(2,3,4)
    plot(EQD.psi_rho_pol, EQD.pressure, 'b', 'linewidth', 2);
    xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Pressure (Newton/m2)')
    
  % Safety factor
  subplot(2,3,5)
    plot(EQD.psi_rho_pol, EQD.q, 'b', 'linewidth', 2);
    xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Safety factor')
    
  % T factor
  subplot(2,3,6)
    plot(EQD.psi_rho_pol, EQD.F, 'b', 'linewidth', 2);
    xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Poloidal current function F (m Tesla)')        
  


% ----------------------------------------------------------------------------------------------------  
% PLOT B summary
% ----------------------------------------------------------------------------------------------------
  
 figure (5, 'position', [100 100 700 700])
  
  subplot(2,2,1)
  
    % Plot the limiter
    hold on
    
    cmin = mean(EQD.B.POL.R(:)) - 2*std(EQD.B.POL.R(:));
    cmax = mean(EQD.B.POL.R(:)) + 2*std(EQD.B.POL.R(:));
     
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R); 
    caxis([cmin cmax]);
    shading flat
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R, 20, 'k')
    title ('B_{\theta,R}(R,Z)')
    h = colorbar;
    set(h, 'title', '(T)')  
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);    
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'w', 'linewidth', 2);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off
    
 
  subplot(2,2,2)
  
    % Plot the limiter
    hold on
    
    cmin = mean(EQD.B.POL.Z(:)) - 1*std(EQD.B.POL.Z(:));
    cmax = mean(EQD.B.POL.Z(:)) + 1*std(EQD.B.POL.Z(:));
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z);
        caxis([cmin cmax]);
    shading flat
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z, 20, 'k')    
    title ('B_{\theta,Z}(R,Z)')
    h = colorbar;
    %set(h, 'title', '(Weber/m^2)')  
    set(h, 'title', '(T)')  
    
     % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);    
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'w', 'linewidth', 2);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight

    hold off
 
 
  subplot(2,2,3)
  
    idx = 4:3:EQD.nr;
    idy = 4:3:EQD.nz;
    quiver(EQD.R_grid(idx, idy), EQD.Z_grid(idx,idy), EQD.B.POL.R(idx,idy), EQD.B.POL.Z(idx,idy), 1)
    hold on
    xlabel('major radius R (m)', 'fontsize', 12)
    ylabel('Z (m)', 'fontsize', 12)
    title ('B_{\theta}(R,Z)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);    
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis tight
    axis equal
    axis tight
    hold off 
 
 
  subplot(2,2,4)
  
    % Plot the limiter
    hold on

    cmin = mean(EQD.B.TOR(:));
    cmax = max(EQD.B.TOR(:));
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.TOR); 
        caxis([cmin cmax]);
    shading flat
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.TOR, 20, 'k'); 
    title ('B_{\phi}(R,Z)')


    h = colorbar;
    set(h, 'title', '(T)')  
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);    
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'w', 'linewidth', 2);    
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off 

    
% ----------------------------------------------------------------------------------------------------  
% PLOT SUMMARY OF RADIAL PROFILES
% ----------------------------------------------------------------------------------------------------    

figure (6, 'position', [100 100 600 600])

    subplot(2,2,1)
    plot(EQD.R, EQD.psi_1D);
        xlabel('R (m)')
        ylabel('stream function \Psi (Weber/rad)')       
    
    subplot(2,2,2)
    plot(EQD.psi_rho_pol, EQD.F, 'linewidth', 2);
    xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_a - \Psi_0|')
    ylabel('Free function F (m Tesla)')       
    
    subplot(2,2,3)
    plot(EQD.psi_rho_pol, EQD.q, 'linewidth', 2);
    xlabel('s = \surd |\Psi - \Psi_0|/|\Psi_a - \Psi_0|')
    ylabel('Safety factor q')
    
    subplot(2,2,4)
        plot(EQD.R, EQD.B.POL.R_radial, 'linewidth', 2, ...
         EQD.R, EQD.B.POL.Z_radial, 'linewidth', 2, ...
         EQD.R, EQD.B.TOR_radial, 'linewidth', 2)
    hl = legend('B_{\theta,R}(R,Z_{0})', 'B_{\theta,Z}(R,Z_{0})', 'B_{\phi}(R,Z_{0})');
    set(hl, 'box', 'off', 'location', 'southeast')
    xlabel('R (m)')
    ylabel('Magnetic Field (T)')
    axis([0 2 -5 1])
 

 
 
% ----------------------------------------------------------------------------------------------------  
% PLOT B divergence
% ---------------------------------------------------------------------------------------------------- 
figure(8)
pcolor(EQD.R_grid, EQD.Z_grid, log10(abs(EQD.DIVB))); 
colormap jet
shading flat
h = xlabel('R (m)');   set(h, 'fontsize', 12);
h = ylabel('Z (m)');  set(h, 'fontsize', 12);
hc = colorbar;
set(hc, 'title', 'log_{10}(|\nabla B(R,Z)|)', 'fontsize', 12);
axis equal
axis tight
hold on
%hc = contour(RR,ZZ, EQLST.poloidalFlux, 20);
%set(hc,'w', 'linewidth', 1); 
hold off 
 




% ----------------------------------------------------------------------------------------------------  
% PLOT B fields
% ----------------------------------------------------------------------------------------------------
  
 figure (2, 'position', [100 100 1400 800])
  
  subplot(2,3,1)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R, 20, 'k')
    title ('Radial Component of B_{pol}')
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off
    
 
  subplot(2,3,2)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z, 20, 'k')    
    title ('Vertical Component of B_{pol}')
    cmin = min(min(EQD.B.POL.Z));
    cmax = max(max(EQD.B.POL.Z(:,5:end)));
    caxis ([cmin cmax])
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight

    hold off
 
 
  subplot(2,3,3)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.T); 
    shading interp
    title ('|B_{pol}|')
        cmin = min(min(EQD.B.POL.T));
    cmax = max(max(EQD.B.POL.T(:,5:end)));
    caxis ([cmin cmax])
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
   
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off 
 
 
  subplot(2,3,4)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.TOR); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.TOR, 20, 'k'); 
    title ('B_{tor}')
    cmin = min(min(EQD.B.TOR(:,100:end)));
    cmax = max(max(EQD.B.TOR));
%      %caxis ([cmin cmax])
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off 
 
 
subplot(2,3,5)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.TOT); 
    shading interp
    title ('|B|')
    cmin = min(min(EQD.B.TOT));
    cmax = max(max(EQD.B.TOT(:,100:end)));
    %caxis ([cmin cmax])    
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off 
 
 subplot(2,3,6)
    plot(EQD.R, EQD.B.POL.R_radial, 'r--', 'linewidth', 2, ...
         EQD.R, EQD.B.POL.Z_radial, 'b--', 'linewidth', 2, ...
         EQD.R, EQD.B.POL.T_radial, 'b', 'linewidth', 2, ...
         EQD.R, EQD.B.TOR_radial, 'm', 'linewidth', 2, ...
         EQD.R, EQD.B.TOT_radial, 'k', 'linewidth',2)
    legend('Bpol_R', 'Bpol_Z', 'Bpol', 'Btor', '|B|')
    xlabel('R (m)')
    ylabel('Magnetic Field (Weber/m^2)')
    title('Magnetic field profile for Z = Z_{axis}')
    axis([0 2 -5 5])
 
% ----------------------------------------------------------------------------------------------------
% PLOT summary
% ----------------------------------------------------------------------------------------------------
 
% Plot the summary: flux surfaces, pressure, q and current

figure(3, 'position', [100 100 800 800])

    % Plot the equilibrium
    subplot(2,2,1)
        hold on
        for k = 1:2:length(EQD.flux_surface_psi)
        plot(EQD.flux_surface_psi(k).R, EQD.flux_surface_psi(k).Z, 'b', 'linewidth', 1.)
        endfor
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);
        
        % Plot the magnetic axis as read from the EQDSK file
        plot(EQD.R_axis, EQD.Z_axis, 'b+', 'linewidth', 2, 'markersize', 20)
        
        % Plot the centroid of 2nd flux surface
        EQD.PSI_R_axis = mean(EQD.flux_surface_psi(1).R);
        EQD.PSI_Z_axis = mean(EQD.flux_surface_psi(1).Z);
        plot(EQD.PSI_R_axis, EQD.PSI_Z_axis, 'rx', 'linewidth', 2, 'markersize', 20)
    
        MAST_Geometry(fillyn = 1);    

        axis equal
        xlabel('R (m)')
        ylabel('Z (m)')
        box on
        
        title(EQD.comments)
        hold off
        
  % Stream function
    subplot(2,2,2)
    plot(EQD.R, EQD.psi_1D);
        xlabel('major radius (m)')
        ylabel('stream function \Psi (Weber/rad)')   
        
  % Plasma pressure
    subplot(2,2,3)  
        plot(EQD.psi_rho_pol, EQD.pressure, 'b', 'linewidth', 2);
        xlabel('\rho')
        ylabel('Pressure (Newton/m2)')
    
  % Safety factor
    subplot(2,2,4)
    plot(EQD.psi_rho_pol, EQD.q', 'b', 'linewidth', 2);
        xlabel('\rho')
        ylabel('Safety factor')        
        
    
% ----------------------------------------------------------------------------------------------------  
% PLOT B fields only
% ----------------------------------------------------------------------------------------------------
  
 figure (4, 'position', [100 100 1400 700])
  
  subplot(1,3,1)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R, 20, 'k')
    title ('Radial Component of B_{pol}')
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off
    
 
  subplot(1,3,2)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z, 20, 'k')    
    title ('Vertical Component of B_{pol}')
    cmin = min(min(EQD.B.POL.Z));
    cmax = max(max(EQD.B.POL.Z(:,5:end)));
    caxis ([cmin cmax])
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight

    hold off
 
 
  subplot(1,3,3)
  
    % Plot the limiter
    plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
    
    % Plot the B field
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.TOR); 
    shading interp
    contour(EQD.R_grid, EQD.Z_grid, EQD.B.TOR, 20, 'k'); 
    title ('B_{tor}')
    cmin = min(min(EQD.B.TOR(:,100:end)));
    cmax = max(max(EQD.B.TOR));
%      %caxis ([cmin cmax])
    h = colorbar;
    set(h, 'title', '(Weber/m^2)')  
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal
    axis tight
    hold off 
 

% ----------------------------------------------------------------------------------------------------  
% PLOT Bp quivers
% ----------------------------------------------------------------------------------------------------
figure(7)
idx = 4:3:EQD.nr;
idy = 4:3:EQD.nz;
quiver(EQD.R_grid(idx, idy), EQD.Z_grid(idx,idy), EQD.B.POL.R(idx,idy), EQD.B.POL.Z(idx,idy), 1)
axis tight
axis equal
xlabel('major radius R (m)', 'fontsize', 12)
ylabel('Z (m)', 'fontsize', 12)
title([filename ' Poloidal Magnetic Field'], 'fontweight', 'bold')


% ----------------------------------------------------------------------------------------------------

figure (10, 'position', [100 100 1400 800])
  
  % 2D geometry
  subplot(2,3,1)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the stream function
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.psi_2D); 
    title('Stream Function \Psi', 'fontsize', 12, 'fontweight', 'bold')
    shading interp
    h = colorbar;
    set(h, 'title', '(Weber/rad)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, 30, 'k', 'linewidth', 1.);
    contour(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal

    hold off

  
  % Plot radial Psi 
  subplot(2,3,2)
    plot(EQD.psi_norm, EQD.psi_r, 'b', 'linewidth', 2, EQD.r_over_a, EQD.psi_vs_r_over_a, 'r', 'linewidth', 2, EQD.psi_rho_pol, EQD.psi_r, 'g', 'linewidth', 2);
    hold on
    %for k = 1:length(EQD.flux_surface)
    %  x(k) =  EQD.flux_surface(k).rho;
    %  y(k) =  EQD.flux_surface(k).V;
    %end
    %plot(x, y, 'b+', 'linewidth', 1.)
    legend(strrep ('EQD.psi_norm', '_', ' '), strrep ('EQD.psi_vs_r_over_a', '_', ' '), strrep ('EQD.psi_rho_pol', '_', ' '),'location', 'southwest')
    hold off
    
    
    xlabel('\Psi_{norm} = |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|, r/a,  \surd {\Psi_{norm} }')
    ylabel('Stream Function \Psi (Weber/rad)')
    title(strrep (filename, '_', ' '))
    
  subplot(2,3,3)
    plot(EQD.R, EQD.psi_1D, 'b', 'linewidth', 2, EQD.psi_1D_R, EQD.psi_1D_levels, 'r+', 'linewidth', 1.5, EQD.R, EQD.psi_edge*ones(size(EQD.R)), 'k');
    xlabel('radius (m)')
    ylabel(['Stream Function \Psi (Weber/rad) at Z = ' num2str(EQD.Z_axis) ' m'])
      
    
  % Plasma pressure
  subplot(2,3,4)
    plot(EQD.psi_rho_pol, EQD.pressure, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Pressure (Newton/m2)')
    
  % Safety factor
  subplot(2,3,5)
    plot(EQD.psi_rho_pol, EQD.q, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Safety factor')
    
  % T factor
  subplot(2,3,6)
    plot(EQD.psi_rho_pol, EQD.F, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Poloidal current function F (m Tesla)')        
    
    
    


% ---------------------------------------------------------------------------------------------------- 

figure (2, 'position', [100 100 1400 800])
  
  % 2D geometry
  subplot(2,3,1)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the poloidal flux
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D); 
    shading interp
    title('Poloidal Flux \Psi_{p}', 'fontsize', 12, 'fontweight', 'bold')
    h = colorbar;
    set(h, 'title', '(Weber)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, 30, 'k', 'linewidth', 1.);
    contour(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal

    hold off

  
  % Plot radial Psi 
  subplot(2,3,2)
    plot(EQD.PSI_norm, EQD.PSI_r, 'b', 'linewidth', 2, EQD.r_over_a, EQD.PSI_vs_r_over_a, 'r', 'linewidth', 2, EQD.PSI_rho_pol, EQD.PSI_r, 'g', 'linewidth', 2);
    hold on
    %for k = 1:length(EQD.flux_surface)
    %  x(k) =  EQD.flux_surface(k).rho;
    %  y(k) =  EQD.flux_surface(k).V;
    %end
    %plot(x, y, 'b+', 'linewidth', 1.)
    legend(strrep ('EQD.PSI_norm', '_', ' '), strrep ('EQD.PSI_vs_r_over_a', '_', ' '), strrep ('EQD.PSI_rho_pol', '_', ' '),'location', 'southwest')
    hold off
    
    
    xlabel('\Psi_{norm} = |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|, r/a,  \surd {\Psi_{norm} }')
    ylabel('Poloidal Flux (Weber/rad)')
    title(strrep (filename, '_', ' '))
    
  subplot(2,3,3)
    plot(EQD.R, EQD.PSI_1D, 'b', 'linewidth', 2, EQD.PSI_1D_R, EQD.PSI_1D_levels, 'r+', 'linewidth', 1.5, EQD.R, EQD.PSI_edge*ones(size(EQD.R)), 'k');
    xlabel('radius (m)')
    ylabel(['Poloidal Flux \Psi_{p} (Weber) at Z = ' num2str(EQD.Z_axis) ' m'])
      
    
  % Plasma pressure
  subplot(2,3,4)
    plot(EQD.PSI_rho_pol, EQD.pressure, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Pressure (Newton/m2)')
    
  % Safety factor
  subplot(2,3,5)
    plot(EQD.PSI_rho_pol, EQD.q, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Safety factor')
    
  % T factor
  subplot(2,3,6)
    plot(EQD.PSI_rho_pol, EQD.F, 'b', 'linewidth', 2);
    xlabel('\Psi_{norm} = \surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
    ylabel('Poloidal current function F (m Tesla)')        
        
