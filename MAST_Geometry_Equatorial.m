function [] = MAST_Geometry_Equatorial()
% Radii
r_ves = 2;
r_mag = 0.9254;
r_lcfs = [0.231, 1.419];
r_fs = [0.815, 1.04, 0.443, 1.329, 0.638, 1.19];
theta = linspace(0,2*pi,360);
x_mag = r_mag*cos(theta);
y_mag = r_mag*sin(theta);
x_lcfs1 = r_lcfs(1)*cos(theta);
y_lcfs1 = r_lcfs(1)*sin(theta);
x_lcfs2 = r_lcfs(2)*cos(theta);
y_lcfs2 = r_lcfs(2)*sin(theta);
x_fs1 = r_fs(1)*cos(theta);
y_fs1 = r_fs(1)*sin(theta);
x_fs2 = r_fs(2)*cos(theta);
y_fs2 = r_fs(2)*sin(theta);
x_fs3 = r_fs(3)*cos(theta);
y_fs3 = r_fs(3)*sin(theta);
x_fs4 = r_fs(4)*cos(theta);
y_fs4 = r_fs(4)*sin(theta);
x_fs5 = r_fs(5)*cos(theta);
y_fs5 = r_fs(5)*sin(theta);
x_fs6 = r_fs(6)*cos(theta);
y_fs6 = r_fs(6)*sin(theta);
x_ves = 2*cos(theta);
y_ves = 2*sin(theta);


% Portholes from CAD drawing
p1 = [-0.646 2.015];
p2 = [-0.630 1.917];
p3 = [-0.045 2.116];
p4 = [-0.028 2.0];
p5 = [-0.02 2.12];
p6 = [-0.671 2.011];

% NBI dumps
scf = 2/700;
BD1 = scf*[182 565];
BD2 = scf*[200 590];
BD3 = scf*[402 537];
BD4 = scf*[396 512];

% beam
NBI1 = 0.95*scf*[133 -720];
NBI2 = 0.92*scf*[228 602];
NBI3 = 0.95*scf*[288 -676];
NBI4 = 0.92*scf*[418 555];



% Generate the second beam dump
bd_angle = pi/3;
BD5(1) = BD1(1)*cos(bd_angle) + BD1(2)*sin(bd_angle);
BD5(2) = BD1(2)*cos(bd_angle) - BD1(1)*sin(bd_angle);
BD6(1) = BD2(1)*cos(bd_angle) + BD2(2)*sin(bd_angle);
BD6(2) = BD2(2)*cos(bd_angle) - BD2(1)*sin(bd_angle);
BD7(1) = BD3(1)*cos(bd_angle) + BD3(2)*sin(bd_angle);
BD7(2) = BD3(2)*cos(bd_angle) - BD3(1)*sin(bd_angle);
BD8(1) = BD4(1)*cos(bd_angle) + BD4(2)*sin(bd_angle);
BD8(2) = BD4(2)*cos(bd_angle) - BD4(1)*sin(bd_angle);

NBI5(1) = NBI1(1)*cos(bd_angle) + NBI1(2)*sin(bd_angle);
NBI5(2) = NBI1(2)*cos(bd_angle) - NBI1(1)*sin(bd_angle);
NBI6(1) = NBI2(1)*cos(bd_angle) + NBI2(2)*sin(bd_angle);
NBI6(2) = NBI2(2)*cos(bd_angle) - NBI2(1)*sin(bd_angle);
NBI7(1) = NBI3(1)*cos(bd_angle) + NBI3(2)*sin(bd_angle);
NBI7(2) = NBI3(2)*cos(bd_angle) - NBI3(1)*sin(bd_angle);
NBI8(1) = NBI4(1)*cos(bd_angle) + NBI4(2)*sin(bd_angle);
NBI8(2) = NBI4(2)*cos(bd_angle) - NBI4(1)*sin(bd_angle);

%figure(3)

hold on
%line([NBI1(1) NBI2(1)], [NBI1(2) NBI2(2)], 'color', 'b', 'linewidth', 2)
%line([NBI3(1) NBI4(1)], [NBI3(2) NBI4(2)], 'color', 'b', 'linewidth', 2)
%line([NBI5(1) NBI6(1)], [NBI5(2) NBI6(2)], 'color', 'b', 'linewidth', 2)
%line([NBI7(1) NBI8(1)], [NBI7(2) NBI8(2)], 'color', 'b', 'linewidth', 2)

%{
line([BD1(1) BD2(1)], [BD1(2) BD2(2)], 'color', 'b', 'linewidth', 2)
line([BD2(1) BD3(1)], [BD2(2) BD3(2)], 'color', 'b', 'linewidth', 2)
line([BD3(1) BD4(1)], [BD3(2) BD4(2)], 'color', 'b', 'linewidth', 2)
line([BD4(1) BD1(1)], [BD4(2) BD1(2)], 'color', 'b', 'linewidth', 2)
line([BD5(1) BD6(1)], [BD5(2) BD6(2)], 'color', 'b', 'linewidth', 2)
line([BD6(1) BD7(1)], [BD6(2) BD7(2)], 'color', 'b', 'linewidth', 2)
line([BD7(1) BD8(1)], [BD7(2) BD8(2)], 'color', 'b', 'linewidth', 2)
line([BD8(1) BD5(1)], [BD8(2) BD5(2)], 'color', 'b', 'linewidth', 2)

fill([BD1(1) BD2(1) BD3(1) BD4(1) BD1(1)], [BD1(2) BD2(2) BD3(2) BD4(2) BD1(2)], [0, 0, 0.5])
fill([BD5(1) BD6(1) BD7(1) BD8(1) BD5(1)], [BD5(2) BD6(2) BD7(2) BD8(2) BD5(2)], [0, 0, 0.5])

line([NBI1(1) NBI2(1) NBI4(1) NBI3(1) NBI1(1)], [NBI1(2) NBI2(2) NBI4(2) NBI3(2) NBI1(2)], 'color', 'b', 'linewidth', 4);
line([NBI5(1) NBI6(1) NBI8(1) NBI7(1) NBI5(1)], [NBI5(2) NBI6(2) NBI8(2) NBI7(2) NBI5(2)], 'color', 'b', 'linewidth', 4);
%h = fill([NBI1(1) NBI2(1) NBI4(1) NBI3(1) NBI1(1)], [NBI1(2) NBI2(2) NBI4(2) NBI3(2) NBI1(2)], [0.8, 0.8, 0.9]);
%set(h,'EdgeColor','None');
%h = fill([NBI5(1) NBI6(1) NBI8(1) NBI7(1) NBI5(1)], [NBI5(2) NBI6(2) NBI8(2) NBI7(2) NBI5(2)], [0.8, 0.8, 0.9]);
%set(h,'EdgeColor','None');
%}

plot(x_lcfs1, y_lcfs1,'k', 'linewidth', 2)
plot(x_lcfs2, y_lcfs2, 'k', 'linewidth', 2)
plot(x_mag, y_mag, 'k', 'linewidth', 1)
plot(x_fs1, y_fs1, 'k--', 'linewidth', 1)
plot(x_fs2, y_fs2, 'k--', 'linewidth', 1)
plot(x_fs3, y_fs3, 'k--', 'linewidth', 1)
plot(x_fs4, y_fs4, 'k--', 'linewidth', 1)
plot(x_fs5, y_fs5, 'k--', 'linewidth', 1)
plot(x_fs6, y_fs6, 'k--', 'linewidth', 1)
plot(x_ves, y_ves, 'k', 'linewidth', 3)

rot_angle = 0.42331;
k = 1;
for alpha = rot_angle:pi/6:rot_angle+2*pi
  P1(k,1) = p1(1)*cos(alpha) + p1(2)*sin(alpha);
  P1(k,2) = p1(2)*cos(alpha) - p1(1)*sin(alpha);
  P2(k,1) = p2(1)*cos(alpha) + p2(2)*sin(alpha);
  P2(k,2) = p2(2)*cos(alpha) - p2(1)*sin(alpha);
  P3(k,1) = p3(1)*cos(alpha) + p3(2)*sin(alpha);
  P3(k,2) = p3(2)*cos(alpha) - p3(1)*sin(alpha);
  P4(k,1) = p4(1)*cos(alpha) + p4(2)*sin(alpha);
  P4(k,2) = p4(2)*cos(alpha) - p4(1)*sin(alpha);
  P5(k,1) = p5(1)*cos(alpha) + p5(2)*sin(alpha);
  P5(k,2) = p5(2)*cos(alpha) - p5(1)*sin(alpha);
  P6(k,1) = p6(1)*cos(alpha) + p6(2)*sin(alpha);
  P6(k,2) = p6(2)*cos(alpha) - p6(1)*sin(alpha);
  
  line([P1(k,1) P2(k,1)], [P1(k,2) P2(k,2)], 'linewidth', 3)
  line([P3(k,1) P4(k,1)], [P3(k,2) P4(k,2)], 'linewidth', 3)
  line([P5(k,1) P6(k,1)], [P5(k,2) P6(k,2)], 'linewidth', 3)

  k = k + 1;
endfor

axis equal




























