function [EQD] = read_eqdsk(filename, plotyn)
% Function that read the EQDSK equilibria and return a stricture containing the equilibirum
% paramaters
%
% INPUT
% filname	full path to the EQDSK file 
% plotyn	keyword: 1 for plotting, 0 for no plot
%
% OUTOUT
% EQD 		structure containing the EFIT equilibria parameters
%
% Example:
% filename = '/home/andrea/Documents/MASTOrbit/29909EQDSK.eqdsk';
% box
% EQD = read_eqdsk(filename, 1); 
% Reading g files from MAST runs: note that in this case the time of the particular equilibrium
% reconstruction is given in the last 5 digits of the filename: 00210 = 0.210 s

% filename = '/home/andrea/ascot5/python/a5py/a5py/preprocessing/29909EQDSK,eqdsk';


% -----------------------------------------------------------------------------------------
% READ EQDSK FILE
% -----------------------------------------------------------------------------------------

% Open the file for reading
fid = fopen(strcat(filename),'r');

EQD.filename = filename;


  % Read the 1st line
  % Note that the comment line might be different for differnt EQDSK files as there is no
  % standard way to write it: this might impact the code between dotted lines
  
 
% ..........................................................................................
  comment_length = 49;					% length in characters of the comments
  EQD.comments = fscanf(fid,'%c', comment_length); 	% comments 
  mah = strsplit(EQD.comments, {"#", " "});
  EQD.pulseNumber = str2num(mah{3});                    % pulse number
  if (strcmp(filename(1), 'g') == 1)
    EQD.selectedtime = str2num(filename(end-4:end))/1000; % time [s]
  else
    EQD.selectedtime = - 1;
  end
 
  % ..........................................................................................
  EQD.id = fscanf(fid,'%i', 1);    	                % indentification character string
  EQD.nr = fscanf(fid,'%i',1);       	                % number of horizontal R grid points
  EQD.nz = fscanf(fid,'%i',1);      	                % number of vertical Z grid points
  %EQD.solver = fscanf(fid,'%c', 1);    	        % solver

  % Read the 2nd line
  EQD.box_width = fscanf(fid,'%e',1); 	                % bow width [m]
  EQD.box_height = fscanf(fid,'%e',1); 	                % box height [m]
  EQD.rcentre = fscanf(fid,'%e',1); 	                % R in meter of vacuum toroidal magnetic field [m]
  EQD.box_left_rad = fscanf(fid,'%e',1); 	        % left position of the R grid [m]
  EQD.zmid = fscanf(fid,'%e',1); 	     	        % Z centre of computational box in meter [m]

  % Read the 3rd line
  EQD.R_axis = fscanf(fid,'%e',1); 	                % magentic axis R coordinate [m]
  EQD.Z_axis = fscanf(fid,'%e',1);	                % magnetic axis Z coordinate [m]                                         
  EQD.simag = fscanf(fid,'%e',1);	                % psi stream function at axis [Weber/rad]
  EQD.sibry = fscanf(fid,'%e',1);	                % psi stream function at edge [Weber/rad]
  EQD.Btor0 = fscanf(fid,'%e',1);	                % Vacuum toroidal magnetic field at rcentre [Tesla]

  % Read the 4th line
  EQD.current = fscanf(fid,'%e',1);	                % total plasma current [Ampere]
  fscanf(fid,'%e',4);			                % skip the next 4 entries

  % Skip the 5th line
  fscanf(fid,'%e',5);	

  % Poloidal current function F = RBtor 
  EQD.F = fscanf(fid,'%e',EQD.nr);	               % [m Tesla]

  % Pressure P
  EQD.pressure = fscanf(fid,'%e',EQD.nr);	       % [Newton/m2] = [Pascals]

  % Poloidal current function prime
  EQD.FF_prime = fscanf(fid,'%e',EQD.nr);              % [m Tesla^2 / (Weber/rad)]
  
  % Pressure gradient P'
  EQD.pressure_prime = fscanf(fid,'%e',EQD.nr);        % [Newton/m^2 / (Weber/rad)]

  % PSI/2 pi 
  EQD.psi = fscanf(fid,'%e',EQD.nr*EQD.nz);	        %[Weber/rad]

  % Safety factor q-profile
  EQD.q = fscanf(fid,'%e',EQD.nr); 

  % Number of points for the LCFS and limiter boundary
  EQD.LCFS_nop = fscanf(fid,'%5i',1);
  EQD.LIMITER_nop = fscanf(fid,'%5i',1);

  % Radial and vertical coordiantes of the LCFS
  tmp = fscanf(fid,'%e',EQD.LCFS_nop*2);
  EQD.LCFS_R = tmp(1:2:end);			       % [m]
  EQD.LCFS_Z = tmp(2:2:end);			       % [m]
  clear tmp

  % Radial and vertical coordinates of the LIMITER (this is not available in MAST)
  tmp= fscanf(fid,'%e',EQD.LIMITER_nop*2);
  EQD.LIMITER_R = tmp(1:2:end);			       % [m]
  EQD.LIMITER_Z = tmp(2:2:end);			       % [m]
  clear tmp

fclose(fid);
% -----------------------------------------------------------------------------------------



% -----------------------------------------------------------------------------------------
% SPATIAL COORDINATES
% -----------------------------------------------------------------------------------------

% Shift the plasma vertically?
shift_equil = 0;
if (shift_equil == 1)
    EQD.Vertical_shift = 0.0;       % [m]
    EQD.Radial_shift = -0.0;        % [m]
    EQD.zmid = EQD.zmid + EQD.Vertical_shift;
    EQD.Z_axis = EQD.Z_axis + EQD.Vertical_shift;
    EQD.R_axis = EQD.R_axis + EQD.Radial_shift;
    EQD.LCFS_Z = EQD.LCFS_Z + EQD.Vertical_shift;
    EQD.LCFS_R = EQD.LCFS_R + EQD.Radial_shift;
    EQD.box_left_rad = EQD.box_left_rad + EQD.Radial_shift;
end

% Generate 1D radial and vertical arrays of spatial coordinates 
EQD.R = linspace(EQD.box_left_rad, EQD.box_left_rad + EQD.box_width, EQD.nr);
EQD.Z = linspace(EQD.zmid-EQD.box_height/2, EQD.zmid+EQD.box_height/2,EQD.nz);


% Generate the 2D grid of radial points from the 1D arrays
[EQD.R_grid, EQD.Z_grid] = meshgrid(EQD.R, EQD.Z);


% -----------------------------------------------------------------------------------------
% COCOS IDENTIFICATION AND CONVERSION
% -----------------------------------------------------------------------------------------

% Convert the equilibrium from COCOS 5 to COCOS 7
[EQD] = cocos_eqdsk(EQD, convert =0 , printyn = 0);

% Check the conversion
cocos_eqdsk(EQD, convert = 0 , printyn = 0);

% -----------------------------------------------------------------------------------------
% 2D STREAM and POLOIDAL FLUX functions
% -----------------------------------------------------------------------------------------

% Calculate the poloidal flux in Weber
EQD.PSI = 2*pi*EQD.psi;	                        %[Weber]

% Reshape stream function and poloidal flux to produce a 2D grid  
EQD.psi_2D = reshape(EQD.psi,EQD.nr,EQD.nz);    % [Weber/rad]
EQD.psi_2D = EQD.psi_2D';
EQD.PSI_2D = reshape(EQD.PSI,EQD.nr,EQD.nz);    % [Weber]
EQD.PSI_2D = EQD.PSI_2D';


% -----------------------------------------------------------------------------------------
% POLOIDAL MAGNETIC FIELD CALCULATION
% -----------------------------------------------------------------------------------------
diffR =  diff(EQD.R);
diffZ =  diff(EQD.Z);
   
% Calculates the magnetic field components from the stream function
[EQD.dpsi_dR, EQD.dpsi_dz] = gradient(EQD.psi_2D, diffR(1), diffZ(1));  % [Weber m⁻1]

if (EQD.COCOS.ID == 3 || EQD.COCOS.ID == 7)
    EQD.B.POL.R = (-1./EQD.R).*EQD.dpsi_dz;	        % [Weber m⁻2]
    EQD.B.POL.Z = (1./EQD.R).*EQD.dpsi_dR;	        % [Weber m⁻2]
elseif (EQD.COCOS.ID == 1 || EQD.COCOS.ID == 5)
    EQD.B.POL.R = (1./EQD.R).*EQD.dpsi_dz;	        % [Weber m⁻2]
    EQD.B.POL.Z = (-1./EQD.R).*EQD.dpsi_dR;	        % [Weber m⁻2]
end  
    
% Total poloidal magnetic field    
EQD.B.POL.T = sqrt(EQD.B.POL.R.^2 + EQD.B.POL.Z.^2);	% [Weber m⁻2]


% -----------------------------------------------------------------------------------------
% STREAM FUNCTION AND POLOIDAL FLUX CALCULATION
% -----------------------------------------------------------------------------------------

% Calculates the stream function on the LCFS
EQD.LCFS_psi = interp2(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, EQD.LCFS_R, EQD.LCFS_Z);    % [Weber/rad]
EQD.LCFS_psi_ave = mean(EQD.LCFS_psi);                                                 % [Weber/rad]


% Calculates the stream function at the position of R_axis and Z_axis
% (this corresponds to  EQD.simag)
EQD.psi_axis = interp2(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, EQD.R_axis, EQD.Z_axis);    % [Weber/rad]
EQD.PSI_axis = 2*pi*EQD.psi_axis;                                                      % [Weber]


% Calulates the polodail flux at the edge as the average
% of the stream function on the LCFS (this corresponds to  EQD.sibry)
EQD.psi_edge = EQD.LCFS_psi_ave;        % [Weber/rad]
EQD.PSI_edge = 2*pi*EQD.psi_edge;       % [Weber]


% Generate the radial profile of the stream function
Z_AXIS = EQD.Z_axis*ones(size(EQD.R));
EQD.psi_1D = interp2 (EQD.R_grid, EQD.Z_grid, EQD.psi_2D, EQD.R, Z_AXIS);    % [Weber/rad]
EQD.PSI_1D = 2*pi*EQD.psi_1D;                                                % [Weber]


% Normalized radial stream function from [0, psi_w = psi_edge]
EQD.psi_1D_norm = abs(EQD.psi_1D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);
EQD.psi_2D_norm = abs(EQD.psi_2D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);




% Determine the radial coordinate where the 1D profile is equal to the 
% the PSI at the edge
%if (EQD.COCOS.sigma.DPSI == 1)
if (EQD.COCOS.ID == 1 || EQD.COCOS.ID == 5)
    psi_index = find(EQD.psi_1D == min(EQD.psi_1D(10:end)) );	% finds the index of R where PSI is min (check the sign of PSI!)
%elseif (EQD.COCOS.sigma.DPSI == -1)
elseif (EQD.COCOS.ID == 3 || EQD.COCOS.ID == 7)
    psi_index = find(EQD.psi_1D == max(EQD.psi_1D(10:end)) );	% finds the index of R where PSI is max (check the sign of PSI!)
end
EQD.R_edge = interp1(EQD.psi_1D(psi_index:end), EQD.R(psi_index:end), EQD.psi_edge);

% Calculates the normalized stream function coordinate
% Note that the number of points is arbitrary.
EQD.psi_n = EQD.nr;
EQD.psi_r = linspace(EQD.psi_axis, EQD.psi_edge, EQD.psi_n);	
EQD.psi_norm = (abs(EQD.psi_r-EQD.psi_axis)./abs(EQD.psi_edge-EQD.psi_axis));
EQD.psi_rho_pol = sqrt(EQD.psi_norm);
EQD.rho2D = sqrt((abs(EQD.psi_2D-EQD.psi_axis)./abs(EQD.psi_edge-EQD.psi_axis)));




% Calculate the radial coordinates corresponding to psi_r
%if (EQD.COCOS.sigma.DPSI == 1)
if (EQD.COCOS.ID == 1 || EQD.COCOS.ID == 5)
    idx = find(EQD.psi_1D <= EQD.psi_axis);
%elseif (EQD.COCOS.sigma.DPSI == -1)
elseif (EQD.COCOS.ID == 3 || EQD.COCOS.ID == 7)
    idx = find(EQD.psi_1D >= EQD.psi_axis);
end
EQD.psi_R_outboard = interp1(EQD.psi_1D(idx:end), EQD.R(idx:end), EQD.psi_r);
EQD.psi_R_inboard = interp1(EQD.psi_1D(4:idx), EQD.R(4:idx), EQD.psi_r);


% Calculate the q profile vs major radius
temp_r = linspace(min(EQD.R), max(EQD.R), 1000);
temp1 = interp1([EQD.R_axis EQD.psi_R_inboard], [EQD.q(1); EQD.q], temp_r);
idx1 = find( isna(temp1) == 1 );
temp1(idx1) = 0;
temp2 = interp1([EQD.R_axis EQD.psi_R_outboard], [EQD.q(1); EQD.q], temp_r);
idx2 = find( isna(temp2) == 1 );
temp2(idx2) = 0;
EQD.q_radial = interp1(temp_r, temp1 + temp2, EQD.R);
clear temp1 temp2;

% Generate the radial profile from R axis to R edge at regular intervals
% and calculates the stream function at these points to be used for the countour plot
%EQD.psi_levels_no = 50;
%EQD.psi_1D_R = linspace(EQD.R_axis, EQD.R_edge, EQD.psi_levels_no);
EQD.psi_1D_R = linspace(EQD.R_axis, EQD.R_edge, EQD.psi_n);
EQD.psi_1D_levels = interp1(EQD.R, EQD.psi_1D, EQD.psi_1D_R);


% Generates the normalized minor radius r/a where a = R_edge - R_axis along the radius
% for Z = Z_axis from the maximum value of PSI to the edge value
EQD.r_over_a = (EQD.psi_1D_R - EQD.R_axis)/(EQD.R_edge - EQD.R_axis);
EQD.psi_vs_r_over_a = EQD.psi_1D_levels;


% Normalized 2D stream function
EQD.psi_2D_min = min(EQD.psi_2D(:));
EQD.psi_2D_max = max(EQD.psi_2D(:));
EQD.psi_2D_norm = (EQD.psi_2D - EQD.psi_axis)/(EQD.psi_edge - EQD.psi_axis);
EQD.psi_2D_norm_2 = (EQD.psi_2D - EQD.psi_2D_max)/(EQD.psi_2D_min - EQD.psi_2D_max);

% Set to zero those values that are < 0 (there should not be none, numerical roudings)
idx = find(EQD.psi_2D_norm < 0);
EQD.psi_2D_norm(idx) = 0;



% -----------------------------------------------------------------------------------------
% 2D CURRENT AND PRESSURE
% -----------------------------------------------------------------------------------------

% Current flux function and Pressure in 2D
% Note that f and P are defined ONLY within the LCFS!!
EQD.F_2D = zeros(size(EQD.psi_2D));
EQD.P_2D = zeros(size(EQD.psi_2D));
for n = 1:EQD.nz
  EQD.F_2D(n,:) = interp1(EQD.psi_norm, EQD.F, EQD.psi_2D_norm(n,:));      % [m Telsa]
  EQD.P_2D(n,:) = interp1(EQD.psi_norm, EQD.pressure, EQD.psi_2D_norm(n,:));       % [Newton/m2]
end

% Set to zero F and P outside the LCFS
idx = find(isna(EQD.F_2D) == 1);
EQD.F_2D(idx) = EQD.Btor0*EQD.rcentre;
idx = find(isna(EQD.P_2D) == 1);
EQD.P_2D(idx) = 0;


% -----------------------------------------------------------------------------------------
% B FIELDS
% -----------------------------------------------------------------------------------------

% Calculates the toroidal B field
EQD.B.TOR = (1./EQD.R).*EQD.F_2D;                        % [Tesla]

% Calculates the total B field
EQD.B.TOT = sqrt(EQD.B.TOR.^2 + EQD.B.POL.R.^2 + EQD.B.POL.Z.^2);                       % [Tesla]

% Calculates the radial profile of the bfields in [Tesla]
EQD.B.TOR_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.B.TOR, EQD.R, Z_AXIS);
EQD.B.POL.R_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.B.POL.R, EQD.R, Z_AXIS);
EQD.B.POL.Z_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.B.POL.Z, EQD.R, Z_AXIS);
EQD.B.POL.T_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.B.POL.T, EQD.R, Z_AXIS);
EQD.B.TOT_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.B.TOT, EQD.R, Z_AXIS);
EQD.B.TOT_axis =  interp2(EQD.R_grid, EQD.Z_grid, EQD.B.TOT, EQD.R_axis, EQD.Z_axis);

% Calculates the divergence of the equilibrium
[EQD.DIVBPR] = gradient(EQD.R_grid.*EQD.B.POL.R, EQD.R)./EQD.R_grid;  
[EQD.DIVBPZ] = gradient(EQD.B.POL.Z, EQD.R, EQD.Z);
EQD.DIVB = EQD.DIVBPR + EQD.DIVBPZ;              


% -----------------------------------------------------------------------------------------
% Calculates the radial profile of the pressure in [Newton/m2]
% -----------------------------------------------------------------------------------------
EQD.P_radial = interp2(EQD.R_grid, EQD.Z_grid, EQD.P_2D, EQD.R, Z_AXIS);


% -----------------------------------------------------------------------------------------
% NORMALIZED POLOIDAL FLUX COORDINATE
% -----------------------------------------------------------------------------------------

% Calculates for all the grid points the corresponding rho and r/a
EQD.rho = interp1(EQD.psi_r, EQD.psi_norm, EQD.psi);
EQD.rho(isnan(EQD.rho)) = 0;
EQD.rho = reshape(EQD.rho,EQD.nr,EQD.nz);
EQD.rho = EQD.rho';

EQD.roa = interp1(EQD.psi_vs_r_over_a, EQD.r_over_a, EQD.psi);
EQD.roa(isnan(EQD.roa)) = 1;
EQD.roa = reshape(EQD.roa, EQD.nr, EQD.nz);
EQD.roa = EQD.roa';

% -----------------------------------------------------------------------------------------


% -----------------------------------------------------------------------------------------
% FLUX SURFACES COORDINATES
% -----------------------------------------------------------------------------------------
[flux_surface] = flux_surfaces_coordinates(EQD, 0);      %
EQD.flux_surface_psi = flux_surface;
clear flux_surface


% -----------------------------------------------------------------------------------------
% SPLINE INTERPOLATION OF psi TO CALCULATE B FIELDS
% -----------------------------------------------------------------------------------------
N = 200;
EQD.INTERP.RI = linspace(EQD.R_grid(1), EQD.R_grid(end), N);
EQD.INTERP.ZI = linspace(EQD.Z_grid(1), EQD.Z_grid(end), N);
[EQD.INTERP.RR, EQD.INTERP.ZZ] = meshgrid(EQD.INTERP.RI,EQD.INTERP.ZI);
EQD.INTERP.psi_2D = interp2(EQD.R_grid, EQD.Z_grid, EQD.psi_2D, EQD.INTERP.RR, EQD.INTERP.ZZ, 'spline');
EQD.INTERP.PSI_2D = interp2(EQD.R_grid, EQD.Z_grid, EQD.PSI_2D, EQD.INTERP.RR, EQD.INTERP.ZZ, 'spline');
EQD.INTERP.psi_2D_norm = interp2(EQD.R_grid, EQD.Z_grid, EQD.psi_2D_norm, EQD.INTERP.RR, EQD.INTERP.ZZ, 'spline');
EQD.INTERP.F_2D = interp2(EQD.R_grid, EQD.Z_grid, EQD.F_2D, EQD.INTERP.RR, EQD.INTERP.ZZ, 'spline');

diffRI = diff(EQD.INTERP.RI);
diffZI = diff(EQD.INTERP.ZI);
[dWR, dWZ] = gradient(EQD.INTERP.psi_2D, diffRI(1), diffZI(1));
EQD.INTERP.BPR = (1./EQD.INTERP.RI).*dWZ;
EQD.INTERP.BPZ = (-1./EQD.INTERP.RI).*dWR;
EQD.INTERP.BP = sqrt(EQD.INTERP.BPR.^2 + EQD.INTERP.BPZ.^2);
EQD.INTERP.BTOR = (1./EQD.INTERP.RI).*EQD.INTERP.F_2D;
EQD.INTERP.B = sqrt(EQD.INTERP.BPR.^2 + EQD.INTERP.BPZ.^2 + EQD.INTERP.BTOR.^2);

ZI_AXIS = EQD.Z_axis*ones(size(EQD.INTERP.RI));
EQD.INTERP.radial.BTOR   = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.BTOR, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.radial.BPR    = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.BPR, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.radial.BPZ    = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.BPZ, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.radial.BP     = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.BP, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.radial.B      = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.B, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.radial.B_axis = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.B, EQD.R_axis, EQD.Z_axis);

EQD.INTERP.psi_1D = interp2 (EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.psi_2D, EQD.INTERP.RI, ZI_AXIS);
EQD.INTERP.PSI_1D = 2*pi*EQD.INTERP.psi_1D;
EQD.INTERP.psi_axis = interp2(EQD.INTERP.RR, EQD.INTERP.ZZ, EQD.INTERP.psi_2D, EQD.R_axis, EQD.Z_axis);
EQD.INTERP.PSI_axis = 2*pi*EQD.INTERP.psi_axis;                                            

EQD.INTERP.psi_1D_norm = abs(EQD.INTERP.psi_1D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);
EQD.INTERP.psi_2D_norm = abs(EQD.INTERP.psi_2D-EQD.psi_axis)/abs(EQD.psi_edge - EQD.psi_axis);


% -----------------------------------------------------------------------------------------
% PLOT THE EQUILIBRIUM
% -----------------------------------------------------------------------------------------
if (plotyn == 1)
 plot_eqdsk(EQD, filename);
end
% Convert for using it in FORCE

%   [EQLST] = EQD2EQLST(EQD);
