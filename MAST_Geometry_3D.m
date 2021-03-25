function [] = MAST_Geometry_3D(EQLST);


[MG] = MOrbit_MAST_Geometry(fillyn = 0, plotyn = 0);
hold on

  Nstep = 20;
%{
  % Draws the last closed flux surface
  [x y z] = revolutionSurface([EQLST.rb(1:10:end)' EQLST.zb(1:10:end)'], linspace(0, pi, Nstep));
  LCFScolor = [0.8 0.8 1.0];
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  %set(h, 'facecolor', 'none')
  set(h, 'facecolor', LCFScolor)
 %} 
  % Draws the polidal field coils
  % coilcolor = [0.914 0.90 0.278];
  coilcolor = [0.9 0.9 0.3];
  divertorcolor = [0.7 0.7 0.7]; 
  
  [x y z] = revolutionSurface(MG.Coils.P2U, linspace(0, 2*pi, 2*Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P4U, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P5U, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)  
    [x y z] = revolutionSurface(MG.Coils.P6U, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P2L, linspace(0, 2*pi, 2*Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P4L, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P5L, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)  
    [x y z] = revolutionSurface(MG.Coils.P6L, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
 
  set(h, 'facecolor', coilcolor)
      [x y z] = revolutionSurface(  MG.Divertor.part1, linspace(0, pi, Nstep));
  h = surf(x,y,z);
  
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', divertorcolor)  
      [x y z] = revolutionSurface(  MG.Divertor.part2, linspace(0, 2*pi, 2*Nstep));
  h = surf(x,y,z);
  
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', divertorcolor)    
      [x y z] = revolutionSurface(  MG.Divertor.part3, linspace(0, 2*pi, 2*Nstep));
  h = surf(x,y,z);
 
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', divertorcolor)  
      [x y z] = revolutionSurface(  MG.Divertor.part4, linspace(0, 2*pi, 2*Nstep));
  h = surf(x,y,z);
  return
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', divertorcolor)      
hold off
axis equal
axis off
box off
hidden
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Z position (m)')
view(308.06, 12.883)

hold off
