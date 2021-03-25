function [X, Y, FIDD] = fastions_dummy(filename)

PA = nc_read(filename, 'A_D_NBI');
EN = nc_read(filename, 'E_D_NBI');
FID = nc_read(filename, 'F_D_NBI');
R2D = nc_read(filename, 'R2D');
Z2D = nc_read(filename, 'Z2D');



[X Y] = meshgrid (linspace(0,200,200), linspace(-150, 150, 200));


% Interpolates the FID on the grid for plotting
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);
w = sum(sum(FID,3),2)*DPA*DEN;
FIDD = griddata(R2D, Z2D, w, X, Y, 'linear');
u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u


