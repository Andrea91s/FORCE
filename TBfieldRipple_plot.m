% Makes plot of the TF ripple


% Grid for the plot
R = linspace(0.23, 2, 200);
P = linspace(0,2*pi, 360);
[RR, PP] = meshgrid(R,P);
xx = RR.*cos(PP);
yy = RR.*sin(PP);
zz = zeros(size(xx));
Re = 1.5;       % plasma edge radius (m)

% 2D TF ripple
Ro = 0.96;      % magnetic axis (m)
Bto = 0.55;     % vacuum toroidal field on axis (T)

% Calculates the 2D TF ripple field
[TFBx, TFBy, TFBz, TFB, TFBphi, TFBR] = TBfieldRipple(xx,yy,zz, Ro, Bto);

figure(1)
MAST_Geometry_Equatorial;
hold on
    pcolor(xx,yy, TFBR); 
    shading interp; 
    plot(Ro*cos(P), Ro*sin(P), 'r--', 'linewidth', 1)
    plot(Re*cos(P), Re*sin(P), 'b--', 'linewidth', 1)
hold off

h = colorbar;
set(h, 'title', 'B_{ripple,R}/B_{\phi}(R_0)', 'fontsize', 12)
xlabel('X (m)', 'fontsize', 12)
ylabel('Y (m)', 'fontsize', 12)
xlim([0 2])
title(['B_{ripple,R}(R,\phi), R_0 = ' num2str(Ro) ' m, B_{\phi}(R_0) = ' num2str(Bto) ' T'], 'fontsize', 12)

return

% Calculates the TF ripple field at different radial locations


phi = linspace(0, 2*pi, 500);
xe = Re*cos(phi);
ye = Re*sin(phi);
ze = zeros(size(xe));
[TFBxe, TFBye, TFBze, TFBe, TFBphie, TFBRe] = TBfieldRipple(xe,ye,ze, Ro, Bto);

xc = Ro*cos(phi);
yc = Ro*sin(phi);
zc = zeros(size(xc));
[TFBxc, TFByc, TFBzc, TFBc, TFBphic, TFBRc] = TBfieldRipple(xc,yc,zc, Ro, Bto);

figure(2)
plot(phi, TFBRe/Bto, 'b', 'linewidth', 2, phi, TFBRc/Bto, 'r', 'linewidth', 2)
xlabel('\phi (rad)', 'fontsize', 12)
ylabel('B_{ripple,R}(R, \phi)/B_{\phi}(R_0)', 'fontsize', 12)
title(['B_{ripple,R}(R,\phi), R_0 = ' num2str(Ro) ' m, B_{\phi}(R_0) = ' num2str(Bto) ' T'], 'fontsize', 12)
legend(['R = ' num2str(Re) ' m'], ['R = ' num2str(Ro) ' m'])
axis([0 2*pi -0.015 0.015])
