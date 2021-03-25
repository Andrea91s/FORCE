function [flux_surface] = flux_surfaces_coordinates(EQD, plotyn)
% Function that generate flux surfaces from the poloidal flux
%
% INPUT
% EQD		strcuture containing the poloidal flux
% plotyn	keyword: 1 for plotting, 0 for no plot
% psi_select    0 to select psi [Weber/rad]
%               1 to select PSI [Weber]                
%
% OUTPUT
% flux_surface  structure with NLEV elements
%               .R radial coordinate [m]
%               .Z vertical coordinate [m]
%               .V poloida flux value [Weber/rad]
%               .Rmax maximum radial position of each flux surface [m]
%               .Zmax vertical coordinate of the maximum radial position of each flux surface [m]
%               .rho normalized poloidal flux



% Choose the right
%switch (psi_select)
%    case 0
%        EQD.psi_1D_levels(end) = EQD.psi_1D_levels(end)*1.000;
%        [C, LEV] = contourc (EQD.R_grid, EQD.Z_grid, EQD.psi_2D, flipdim(EQD.psi_1D_levels,2));
%        psi_min_index = find(EQD.psi_1D == min(EQD.psi_1D) );	% finds the index of R where PSI is min
%    case 1
%        EQD.PSI_1D_levels(end) = EQD.PSI_1D_levels(end)*1.000;
%        [C, LEV] = contourc (EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, flipdim(EQD.PSI_1D_levels,2));  
%        PSI_min_index = find(EQD.PSI_1D == min(EQD.PSI_1D) );	% finds the index of R where PSI is min
%endswitch

EQD.psi_1D_levels(end) = EQD.psi_1D_levels(end)*1.000;
%[C, LEV] = contourc (EQD.R_grid, EQD.Z_grid, EQD.psi_2D, flipdim(EQD.psi_1D_levels,2));
%psi_min_index = find(EQD.psi_1D == min(EQD.psi_1D) );	% finds the index of R where PSI is min

[C, LEV] = contourc (EQD.R_grid, EQD.Z_grid, EQD.psi_2D, flip(EQD.psi_1D_levels,2));

%if (EQD.COCOS.sigma.DPSI == 1)
if (EQD.COCOS.ID == 1 || EQD.COCOS.ID == 5)
    psi_index = find(EQD.psi_1D == min(EQD.psi_1D) );	% finds the index of R where PSI is min
%elseif (EQD.COCOS.sigma.DPSI == -1)
elseif (EQD.COCOS.ID == 3 || EQD.COCOS.ID == 7)
    psi_index = find(EQD.psi_1D == max(EQD.psi_1D) );	% finds the index of R where PSI is max
endif

% Keep the flux at the LCFS slightly above its nominal value to avoid
% the edge flux surface to be computed open! This is necessary for 
% scenario four!
% EQD.psi_1D_levels(end) = EQD.psi_1D_levels(end)*1.000;

% Compute the spatial coordinates of the contour lines of PSI on the levels specified in the
% EQD structure
% [C, LEV] = contourc (EQD.R_grid, EQD.Z_grid, EQD.psi_2D, flipdim(EQD.psi_1D_levels,2));



% Build a structure containing the x,y coordinates of the equilibrium
l = 1;
o = 1;
while (l < length(C(2,:)))
  m(o) = C(2,l);	% number of elements in contour level o
  n = l+1;		% 1st index of contour level o
  curve(o).x = C(1, n:n+m(o)-1);
  curve(o).y = C(2, n:n+m(o)-1);
  curve(o).value = C(1,n-1);
  l = n + m(o);
  o = o + 1;
endwhile
O = o - 1; % number of levels generated
%fprintf('Number of level generated: %d\n', O);

% remove the non-closed open contour lines
DT = 1E-10;     % difference threshold: if DX, DY > DT the flux surface is not closed
m = 1;
for k = 1:O
  %printf('\n\nIso-level %d at Psi = %f\n', k, curve(k).value); fflush(stdout());
  %printf('x_start = %f \tx_end = %f\ty_start = %f\ty_end = %f\n', curve(k).x(1), curve(k).x(end),curve(k).y(1), curve(k).y(end)); fflush(stdout());
  DX = abs(curve(k).x(1) - curve(k).x(end));
  DY = abs(curve(k).y(1) - curve(k).y(end));
  if (DX < DT && DY< DT)
    %printf('\tClosed iso-level no %d at Psi = %f\n', m, curve(k).value); fflush(stdout());
    %printf('\tx_start = %f \tx_end = %f\ty_start = %f\ty_end = %f\n', curve(k).x(1), curve(k).x(end),curve(k).y(1), curve(k).y(end)); fflush(stdout());
    flux_surface(m).R = curve(k).x;
    flux_surface(m).Z = curve(k).y;
    flux_surface(m).V = curve(k).value;
    flux_surface(m).Rmax = max(flux_surface(m).R);
    flux_surface(m).Zmax = flux_surface(m).Z (find( flux_surface(m).R ==  flux_surface(m).Rmax));
    %fprintf('Level # %d\nCLOSED level generated: %d\n', k, m);
    m = m + 1;
  %elseif
   %fprintf('Open Level # %d:\n', k);
   %fprintf('x(1) = %1.8f x(end) = %1.8f DX = %1.8g\n', curve(k).x(1), curve(k).x(end), DX);
   %fprintf('y(1) = %1.8f y(end) = %1.8f DY = %1.8g\n', curve(k).y(1), curve(k).y(end), DY);
  endif
endfor

% Calculates the coordinates yet again
M = m-1; % number of closed levels generated
%fprintf('Number of CLOSED level generated: %d\n', M);
%PSI_max_index = find(EQD.psi_1D == max(EQD.psi_1D) );	% finds the index of R where PSI is max 
%PSI_min_index = find(EQD.psi_1D == min(EQD.psi_1D) );	% finds the index of R where PSI is min

for k = 1:M;
  flux_surface(k).rho = sqrt(abs(flux_surface(k).V - flux_surface(M).V)/abs(flux_surface(M).V-flux_surface(1).V));
  flux_surface(k).r = abs(flux_surface(k).Rmax - flux_surface(M).Rmax)/abs(flux_surface(M).Rmax-flux_surface(1).Rmax);
  flux_surface(k).r_over_a = abs(flux_surface(k).Rmax - EQD.R_axis)/abs(EQD.R_edge - EQD.R_axis);
  flux_surface(k).r_interpolated = interp1(EQD.psi_1D(psi_index:end), EQD.R(psi_index:end), flux_surface(k).V);
end


%switch (psi_select)
%    case 0
%        for k = 1:M
%        %flux_surface(k).r_interpolated = interp1(EQD.psi_1D(PSI_max_index:end), EQD.R(PSI_max_index:end), flux_surface(k).V);
%        flux_surface(k).r_interpolated = interp1(EQD.psi_1D(psi_min_index:end), EQD.R(psi_min_index:end), flux_surface(k).V);
%        end
%    case 1
%        for k = 1:M
%        %flux_surface(k).r_interpolated = interp1(EQD.psi_1D(PSI_max_index:end), EQD.R(PSI_max_index:end), flux_surface(k).V);
%        flux_surface(k).r_interpolated = interp1(EQD.PSI_1D(PSI_min_index:end), EQD.R(PSI_min_index:end), flux_surface(k).V);
%        end
%endswitch        




% remove the closed flux surfaces not centred on the magnetic axis
mag_axis_r = mean(flux_surface(1).R);
mag_axis_z = mean(flux_surface(1).Z);
p = 1;
for k = 1:m-1
    if (inpolygon(mag_axis_r, mag_axis_z, flux_surface(k).R, flux_surface(k).Z) == 1)
      flux_surfaces(p).R = flux_surface(k).R;
      flux_surfaces(p).Z = flux_surface(k).Z;
      flux_surfaces(p).V = flux_surface(k).V;
      flux_surfaces(p).rho = flux_surface(k).rho;
      flux_surfaces(p).r = flux_surface(k).r;
      flux_surfaces(p).r_over_a = flux_surface(k).r_over_a;
      flux_surfaces(p).r_interpolated = flux_surface(k).r_interpolated;
      p = p + 1;
    endif
endfor

clear flux_surface;

flux_surface = flux_surfaces;
M = p - 1;



% Plot the results
if (plotyn == 1)
figure(20)
  hold on
    %for k = 1:O
    %  plot(curve(k).x, curve(k).y, 'r', 'linewidth', 1.)
    %endfor
    for k = 1:2:M
      plot(flux_surface(k).R, flux_surface(k).Z, 'b', 'linewidth', 1.)
    endfor
    %for k = 1:p-1
    %  plot(flux_surfaces(k).R, flux_surfaces(k).Z, 'b', 'linewidth', 1.)
    %endfor
    
    % Plot the limiter
    % plot(EQD.LIMITER_R, EQD.LIMITER_Z, 'k', 'linewidth', 2);
    
    % Plot the last closed flux surface
    plot(EQD.LCFS_R, EQD.LCFS_Z, 'r', 'linewidth', 2);
    
    % Plot the magnetic axis as read from the EQDSK file
    plot(EQD.R_axis, EQD.Z_axis, 'b+', 'linewidth', 2, 'markersize', 20)
    
    % Plot the centroid of 2nd flux surface
    EQD.PSI_R_axis = mean(flux_surface(end-1).R);
    EQD.PSI_Z_axis = mean(flux_surface(end-1).Z);
    plot(EQD.PSI_R_axis, EQD.PSI_Z_axis, 'rx', 'linewidth', 2, 'markersize', 20)
  
    MAST_Geometry(fillyn = 1);    

    axis equal
    xlabel('R (m)')
    ylabel('Z (m)')
    box on
  hold off
  
%figure(21)
%  for k = 1:M
%    x(k) =  flux_surface(k).rho;
%    y(k) =  flux_surface(k).V;
%  end
%  plot(EQD.psi_rho, EQD.psi_r, 'ro', 'linewidth', 0.5, x, y, 'b', 'linewidth', 2.)
  %xlabel('Radial coordinate (m)')
  %ylabel('\surd |\Psi - \Psi_0|/|\Psi_1 - \Psi_0|')
endif
