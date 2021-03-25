function [EQDPP] = pp_eqdsk(EQD, plotyn)
% Function that interpolates PSI using Octave piece-wise polynomial splines and calculates
% the magnetic fields taking the derivatives of the splined PSI thus avoiding the use of
% the gradient function.
%
% This routine is used only used to prepare the B field and PSI needed by the full orbit code
% and therefore computes only those equil parameters that are needed by it.
%
% INPUT
% filname	EQD strucutured read from read_eqdsk 
% plotyn	keyword: 1 for plotting, 0 for no plot
%
% OUTOUT
% EQDPP.	structure containing the EFIT equilibria parameters from piece-wise poly splines
%
% Example:
% directory = '/home/marco/Documents/MAST/TRANSP/Kick_Model/EFIT/';
% filename = 'g029976.00160';
% EQD = read_eqdsk(directory, filename, 1);
% EQDPP = pp_eqdsk(EQD, plotyn)


% In MAST, PSI varies significantly in the first 4 radial points (R <= 0.15 m) and this
% causes pp-splines to incorrect model PSI for R > 0.1 5m; since these region is not that
% important, EQDPP equilibrium is calculated starting from a spatial grid which is:
%
% (EQD.nr - 4) x EQD.nz
R      = EQD.R(4:EQD.nr);
Z      = EQD.Z;
psi2D  = EQD.psi_2D(:,4:EQD.nr);
F2D    = EQD.F_2D(:,4:EQD.nr);


% Generate a finer spatial grid on which PSI will be pp-splined and the B fields calculated
NP = 500;
RI = linspace(EQD.R(4), EQD.R(EQD.nr), NP);
ZI = linspace(EQD.Z(1), EQD.Z(EQD.nz), NP);
[RR, ZZ] = meshgrid(RI,ZI);

%------------------------------------------------------------------
% PP_SPLINE of PSI
%------------------------------------------------------------------
% Use PP to interpolate PSI along the radial direction
M = numel(Z);
for m = 1:M
    pp = spline(R, psi2D(m,:));
    W(m,:) = ppval(pp, RI);
    clear pp
end

% Use PP to interpolate PSIR along the vertical direction
N = numel(RI);
for n = 1:N
    pp = spline(Z, W(:,n));
    psi2Dpp(:,n) = ppval(pp, ZI);
    clear pp
end

clear W;

%------------------------------------------------------------------
% PP_SPLINE of F
%------------------------------------------------------------------
% Use PP to interpolate F along the radial direction
M = numel(Z);
for m = 1:M
    pp = spline(R, F2D(m,:));
    W(m,:) = ppval(pp, RI);
    clear pp
end

% Use PP to interpolate F along the vertical direction
N = numel(RI);
for n = 1:N
    pp = spline(Z, W(:,n));
    F2Dpp(:,n) = ppval(pp, ZI);
    clear pp
end

clear W;

%------------------------------------------------------------------
% PP_SPLINE of dPSI/dR and dPSI/dZ
%------------------------------------------------------------------
% Use PPD to calculate dPSI/dR along the radial direction
M = numel(ZI);
for m = 1:M
    pp = spline(RI, psi2Dpp(m,:));
    ppd = ppder(pp,1);
    dWR(m,:) = ppval(ppd, RI);    
    clear ppd
end

% Use PPD to calculate dPSI/dZ along the vertical direction
N = numel(RI);
for n = 1:N
    pp = spline(ZI, psi2Dpp(:,n));
    ppd = ppder(pp,1);    
    dWZ(:,n) = ppval(ppd, ZI);
    clear pp
end

%------------------------------------------------------------------
% Copy partial results into EQDPP structure 
%------------------------------------------------------------------
EQDPP.psi_2D = psi2Dpp;
EQDPP.PSI_2D = 2*pi*psi2Dpp;
EQDPP.F_2D = F2Dpp;
EQDPP.R = RI;
EQDPP.Z = ZI;
EQDPP.R_grid = RR;
EQDPP.Z_grid = ZZ;


%------------------------------------------------------------------
% Computes the B field using pp-splines
%------------------------------------------------------------------
EQDPP.B.POL.R   = (1./RI).*dWZ;	                % [Weber m⁻2]
EQDPP.B.POL.Z   = (-1./RI).*dWR;	        % [Weber m⁻2]
EQDPP.B.TOR     = (1./EQDPP.R).*EQDPP.F_2D;     % [Tesla]

% Calculates the total B field

EQDPP.B.TOT = sqrt(EQDPP.B.TOR.^2 + EQDPP.B.POL.R.^2 + EQDPP.B.POL.Z.^2); 

% but computes the divergence numerically (this is not needed
% by the full orbit code)
[EQDPP.DIVBPR] = gradient(RR.*EQDPP.B.POL.R, RI)./RR;  
[EQDPP.DIVBPZ] = gradient(EQDPP.B.POL.Z,RI,ZI);
EQDPP.DIVB = EQDPP.DIVBPR + EQDPP.DIVBPZ;               % div B



%------------------------------------------------------------------
% Computes the B field using the gradient function intead
%------------------------------------------------------------------
DR = RI(2) - RI(1);
DZ = ZI(2) - ZI(1);
[dpsi2Dpp_dR, dpsi2Dpp_dZ] = gradient(psi2Dpp, DR, DZ);
BPR = (1./RI).*dpsi2Dpp_dZ;	        % [Weber m⁻2]
BPZ = (-1./RI).*dpsi2Dpp_dR;	        % [Weber m⁻2]

DIVBPR = gradient(RR.*BPR, RI)./RR;  
[dummy, DIVBPZ] = gradient(BPZ,RI,ZI);
DIVB = DIVBPR + DIVBPZ;               % div B



%------------------------------------------------------------------
% Computes the other quantities needed by the full orbit code
%------------------------------------------------------------------

% Generate the radial profile of the stream function
Z_AXIS = EQD.Z_axis*ones(size(RI));
EQDPP.psi_1D = interp2 (RR, ZZ, psi2Dpp, RI, Z_AXIS, 'spline');    % [Weber/rad]
EQDPP.PSI_1D = 2*pi*EQDPP.psi_1D;                                  % [Weber]

% Normalized radial stream function from [0, psi_w = psi_edge]
EQDPP.psi_2D_norm = abs(EQDPP.psi_2D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);
EQDPP.psi_1D_norm = abs(EQDPP.psi_1D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);




% -----------------------------------------------------------------------------------------
% PLOT THE ORIGINAL and PP-SPLINED QUANTITIES side by side for comparison
% -----------------------------------------------------------------------------------------
if (plotyn == 1)

 figure (1)
  
  % Plot the stream function PSI
  subplot(1,2,1)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    cmin = min(EQD.psi_2D(:));
    cmax = max(EQD.psi_2D(:));
        
    % Plot the stream function
    pcolor(EQD.R_grid, EQD.Z_grid, EQD.psi_2D); 
    title('Stream Function \Psi', 'fontsize', 12, 'fontweight', 'bold')
    caxis ([cmin cmax])
    colormap jet
    shading flat
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

  subplot(1,2,2)
    % Plot the limiter
    %plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    hold on
    
    % Plot the stream function
    pcolor(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.psi_2D); 
    caxis ([cmin cmax])
    colormap jet
    title('PP - Stream Function \Psi', 'fontsize', 12, 'fontweight', 'bold')
    shading flat
    h = colorbar;
    set(h, 'title', '(Weber/rad)')

    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r--', 'linewidth', 2);
    
    contour(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.psi_2D, 30, 'k', 'linewidth', 1.);
    contour(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.psi_2D, [EQD.sibry EQD.sibry], 'k--', 'linewidth', 2.);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    axis equal

    hold off    
    
    
    
 figure (2)
    
    subplot(1,2,1)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);

        cmin = min(EQD.B.POL.R(:));
        cmax = max(EQD.B.POL.R(:));
    
        % Plot the B field
        pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R);
        caxis ([cmin cmax])
        colormap jet
        shading flat
        contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R, 20, 'k')
        title ('Radial Component of B_{pol}')
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight
        hold off
   
    subplot(1,2,2)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
        
        % Plot the B field
        pcolor(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.POL.R);
        caxis ([cmin cmax])    
        colormap jet
        shading flat
        contour(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.POL.R, 20, 'k')
        title ('Radial Component of PP B_{pol}')
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight
        hold off
        
    
 figure (3)
    
    subplot(1,2,1)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
        
        % Plot the B field
        pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z); 
        colormap jet
        shading flat
        contour(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z, 20, 'k')    
        title ('Vertical Component of B_{pol}')
        cmin = min(min(EQD.B.POL.Z(:,10:end)));
        cmax = max(EQD.B.POL.Z(:));
        caxis ([cmin cmax])
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight

        hold off
    
    
    subplot(1,2,2)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
        
        % Plot the B field
        pcolor(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.POL.Z); 
        colormap jet
        shading flat
        contour(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.POL.Z, 20, 'k')    
        title ('Vertical Component of PP B_{pol}')
        caxis ([cmin cmax])
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight

        hold off
    
    
    
 figure (4)
    
    subplot(1,2,1)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);

        cmin = 0.1*min(EQD.B.TOR(:));
        cmax = max(EQD.B.TOR(:));        
        
        % Plot the B field
        pcolor(EQD.R_grid, EQD.Z_grid, EQD.B.TOR); 
        caxis ([cmin cmax]);
        colormap jet
        shading flat
        contour(EQD.R_grid, EQD.Z_grid, EQD.B.TOR, 20, 'k'); 
        title ('B_{tor}')
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight
        hold off     
        

        

    subplot(1,2,2)
    
        % Plot the limiter
        plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
        hold on
        
        % Plot the last closed flux surface
        plot(EQD.LCFS_R, EQD.LCFS_Z, 'k--', 'linewidth', 2);
        
        % Plot the B field
        pcolor(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.TOR);
        caxis ([cmin cmax]);        
        colormap jet
        shading flat
        contour(EQDPP.R_grid, EQDPP.Z_grid, EQDPP.B.TOR, 20, 'k'); 
        title ('PP B_{tor}')
        h = colorbar;
        set(h, 'title', '(Weber/m^2)')  
        
        xlabel('R (m)')
        ylabel('Z (m)')
        axis equal
        axis tight
        hold off             
        

 figure(5)
 
    subplot(1,2,1)
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

 
    subplot(1,2,2)
    pcolor(EQDPP.R_grid, EQDPP.Z_grid, log10(abs(EQDPP.DIVB))); 
    colormap jet
    shading flat
    h = xlabel('R (m)');   set(h, 'fontsize', 12);
    h = ylabel('Z (m)');  set(h, 'fontsize', 12);
    hc = colorbar;
    set(hc, 'title', 'PP log_{10}(|\nabla B(R,Z)|)', 'fontsize', 12);
    axis equal
    axis tight
    hold on
    %hc = contour(RR,ZZ, EQLST.poloidalFlux, 20);
    %set(hc,'w', 'linewidth', 1); 
    hold off     
    
endif























