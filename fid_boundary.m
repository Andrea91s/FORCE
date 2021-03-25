function [] = fid_boundary(PA, EN, R2D, Z2D, FID, id, IE, CM, Radial_Position)
% Function used to plot the FID at a particular time and position (R,Z)
% together with the particle orbit topology boundaries
%
% NB    The orbit boundaries are added manually in this scritp!
%
% INPUT
% R2D		array		radial position at which the FID has been calculated in cm
% Z2D		array		vertical position at which the FID has been calculated in cm
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% id		integer		array index of FID which is closest to the requested point in the poloidal plane
% CM 		real		Z max of the FID: 0 for auto
% IE		real		maximum injection energy for plotting
% Radial_Position               R in cm to select the correct boundaries

% Example
% filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U16/29880U16_fi_1.cdf';
% [R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 90, RZ = 0, CM = 0, IE = 6E4, plotyn = 1, saveyn = 0);
% fid_boundary(PA, EN, R2D, Z2D, FID, id, IE, CM);
% fid_boundary(U16pre.PA, U16pre.EN, U16pre.R2D, U16pre.Z2D, U16pre.FID, U16pre.id, U16pre.IE, U16pre.CM);

fprintf('Boundaries for R = %f cm\n', Radial_Position);

% Particle orbit Topology Boundaries 
switch (Radial_Position)
    case 90

        % ----------------------------------------------------------------------------------------------------
        % Orbit Topology Boundaries for pulse 29880 at 0.258 s and 0.263 s, R = 0.9 m Z = 0.0 m
        % ----------------------------------------------------------------------------------------------------
        %{
        B(1).x = [0.25 0.25 0.35 0.45 0.55 0.65];
        B(1).y = [0 7.5 12.5 32.5 57.5 62.5]*1E3;
        B(2).x = [0.05 0.05];
        B(2).y = [0 62.5]*1E3;
        %}
        B(1).x = [-1.05 -0.35 -0.35 -0.25 -0.05 -0.05]
        B(1).y = [2.5 7.5 7.5 12.5 12.5 62.5]*1E3;
        B(2).x = [-0.05 -0.05 0.05 0.05 -0.05]
        B(2).y = [2.5 7.5 7.5 2.5 2.5]*1E3;
        fprintf('Number of boundaries: %d\n', numel(B));

    case 100

        % ----------------------------------------------------------------------------------------------------
        % Orbit Topology Boundaries for 29880 at 0.255 s, R = 1.0, Z = 0.0 m
        % ----------------------------------------------------------------------------------------------------
        B(1).x = [-1.05 -0.85 -0.85 -0.75 -0.75 -0.65 -0.65 -0.55 -0.45 -0.45 -0.35 -0.35];
        B(1).y = [12.50 12.50 17.50 17.50 22.50 22.50 17.50 12.50 12.50  7.50  7.50  0.00]*1E3;
        B(2).x = [0.05 0.05];
        B(2).y = [0.00 62.5]*1E3;
        B(3).x = [0.55 0.55 0.65 0.65 0.55 0.55 0.65 0.65 0.75 0.75];
        B(3).y = [0.00 22.5 22.5 27.5 27.5 32.5 32.5 57.5 57.5 62.5]*1E3;
        fprintf('Number of boundaries: %d\n', numel(B));
        
    case 110

        % ----------------------------------------------------------------------------------------------------
        % Orbit Topology Boundaries for pulse 29880 at 0.258 s, R = 1.1 m Z = 0.0 m
        % ----------------------------------------------------------------------------------------------------
        B(1).x = [-0.45 -0.45 -0.35 -0.35];
        B(1).y = [0 12.5 12.5 62.5]*1E3;
        B(2).x = [-0.45 -0.45 -0.35 -0.35 -0.25 -0.25 -0.15 -0.15 -0.05 -0.05];
        B(2).y = [0 12.5 12.5 17.5 17.5 22.5 22.5 7.5 7.5 0]*1E3;
        B(3).x = [-0.35 -0.35 -0.25 -0.25 -0.15 -0.15 -0.05 -0.05];
        B(3).y = [62.5 17.5 17.5 22.5 22.5 7.5 7.5 0]*1E3;
        B(4).x = [0.05 0.05];
        B(4).y = [0 62.5]*1E3;
        B(5).x = [0.65 0.65 0.75 0.75];
        B(5).y = [0 22.5 22.5 62.5]*1E3;
        fprintf('Number of boundaries: %d\n', numel(B));

    case 120
        
        % ----------------------------------------------------------------------------------------------------
        % Boundaries for pulse 29880 at 0.255 and 0.263 s, R = 1.2 m Z = 0.0 m
        % ----------------------------------------------------------------------------------------------------
        B(1).x = [-0.65 -0.65 -0.55 -0.55 -0.45 -0.45 -0.35 -0.35];
        B(1).y = [ 0.00  7.50  7.50 22.50 22.50 57.50 57.50 62.50]*1E3;
        B(2).x = [-0.55 -0.55 -0.45 -0.45 -0.55];
        B(2).y = [17.50 22.50 22.50 17.50 17.50]*1E3;
        B(3).x = [-0.45 -0.45 -0.35 -0.35 -0.45];
        B(3).y = [42.50 52.50 52.50 42.50 42.50]*1E3;
        B(4).x = [-0.15 -0.15 -0.05 -0.05];
        B(4).y = [62.50 47.50 47.50  0.00]*1E3;
        B(5).x = [0.05  0.05];
        B(5).y = [0.00 62.50]*1E3;
        B(6).x = [0.65 0.65 0.75 0.75 0.85 0.85 0.95 0.95];
        B(6).y = [62.50 52.50 52.50 32.50 32.50 57.50 57.50 62.50]*1E3;
        B(7).x = [0.75 0.75 0.85 0.85 0.95 0.95];
        B(7).y = [0.00 27.50 27.50 57.50 57.50 62.50]*1E3;
        fprintf('Number of boundaries: %d\n', numel(B));        

endswitch





% ----------------------------------------------------------------------------------------------------
% Makes the plot
% ----------------------------------------------------------------------------------------------------
[pa, en] = meshgrid(PA, EN);

figure(2)
    pcolor(pa, en, squeeze(FID(id,:,:))');
    colormap jet
    shading flat
    h = xlabel('pitch'); set(h, 'fontsize', 12);
    h = ylabel('energy eV'); set(h, 'fontsize', 12);
    if (CM != 0)
    axis([-1 1 0 IE 0 CM])
    caxis([0 CM])
    else  
        axis([-1 1 0 IE])
    endif
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
    h = colorbar; set(h, 'title', 'cm^{-3} eV^{-1}');  
    
    % Draw the boundaries
    NB = numel(B);
    hold on
        for nb = 1:NB
            h = stairs(B(nb).x, B(nb).y);
            set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        end
    hold off
    

    
figure(3)
    pcolor(pa, en, log10(squeeze(FID(id,:,:))'));
    shading flat
    xlabel('pitch')
    ylabel('energy eV')
    axis([-1 1 0 6E4])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
    colorbar  
    
    % Draw the boundaries
    % Draw the boundaries
    NB = numel(B);
    hold on
        for nb = 1:NB
            h = stairs(B(nb).x, B(nb).y);
            set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        end
    hold off
% ----------------------------------------------------------------------------------------------------



return


% ----------------------------------------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.263 s, R = 1.1 m Z = 0.0 m
% ----------------------------------------------------------------------------------------------------
B(1).x = [-0.45 -0.45 -0.35 -0.35];
B(1).y = [0 7.5 7.5 62.5];
B(2).x = [-0.45 -0.45 -0.35 -0.35 -0.15 -0.15 -0.05 -0.05];
B(2).y = [0 7.5 7.5 17.5 17.5 7.5 7.5 0];
B(3).x = [-0.35 -0.35 -0.15 -0.15 -0.05 -0.05];
B(3).y = [62.5 17.5 17.5 7.5 7.5 0];
B(4).x = [0.05 0.05];
B(4).y = [0 62.5];
B(5).x = [0.65 0.65 0.75 0.75];
B(5).y = [0 22.5 22.5 62.5];     


B(1).x = [-1.05 -0.85 -0.85 -0.75 -0.75 -0.65 -0.65 -0.55 -0.45 -0.45 -0.35 -0.35];
B(1).y = [12.50 12.50 17.50 17.50 22.50 22.50 17.50 12.50 12.50  7.50  7.50  0.00]*1E3;
B(2).x = [0.05 0.05];
B(2).y = [0.00 62.5]*1E3;
B(3).x = [0.55 0.55 0.65 0.65 0.55 0.55 0.65 0.65 0.75 0.75]+DL;
B(3).y = [0.00 22.5 22.5 27.5 27.5 32.5 32.5 57.5 57.5 62.5]*1E3;     





% ----------------------------------------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.255 and 0.263 s, R = 1.2 m Z = 0.0 m with TF ripples
% ----------------------------------------------------------------------------------------------------
B(1).x = [-0.65 -0.65 -0.55 -0.55 -0.45 -0.45 -0.35 -0.35];
B(1).y = [ 0.00  7.50  7.50 27.50 27.50 57.50 57.50 62.50]*1E3;
B(2).x = [-0.55 -0.55 -0.45 -0.45 -0.55];
B(2).y = [17.50 27.50 27.50 17.50 17.50]*1E3;
B(3).x = [-0.45 -0.45 -0.35 -0.35 -0.45];
B(3).y = [42.50 57.50 57.50 42.50 42.50]*1E3;
B(4).x = [-0.15 -0.15 -0.05 -0.05];
B(4).y = [62.50 57.50 57.50  0.00]*1E3;
B(5).x = [0.05  0.05];
B(5).y = [0.00 62.50]*1E3;
B(6).x = [0.55   0.55  0.65  0.65  0.75  0.75  0.85  0.85  0.95  0.95];
B(6).y = [62.50 52.50 52.50 47.50 47.50 32.50 32.50 57.50 57.50 62.50]*1E3;
B(7).x = [0.75  0.75  0.85  0.85  0.95  0.95];
B(7).y = [0.00 27.50 27.50 57.50 57.50 62.50]*1E3;
