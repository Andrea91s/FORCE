function [] = MOrbit_plot_vtor(particle, idx)
% Function that plot the velocity components at the given time
% This is used to check that I have calculated the toroidal velocity correctly.

figure(200)
hold on
drawCircle(0,0,particle.orbit.position.R(idx),'k')

% Draw the vxy vector
drawVector3d([particle.orbit.position.x(idx) particle.orbit.position.y(idx) 0 ], [particle.orbit.velocity.vx(idx) particle.orbit.velocity.vy(idx) 0]/particle.orbit.velocity.vxy(idx), 'b', 'linewidth', 2); 
drawVector3d([particle.orbit.position.x(idx) particle.orbit.position.y(idx) 0 ], [particle.orbit.velocity.vx(idx) 0 0]/particle.orbit.velocity.vxy(1), 'b', 'linewidth', 1); 
drawVector3d([particle.orbit.position.x(idx) particle.orbit.position.y(idx) 0 ], [0 particle.orbit.velocity.vy(idx) 0]/particle.orbit.velocity.vxy(1), 'b', 'linewidth', 1); 

% Draw the toroidal velocity
drawVector3d([particle.orbit.position.x(idx) particle.orbit.position.y(idx) 0 ], [particle.orbit.velocity.vtx(idx) particle.orbit.velocity.vty(idx) 0]/particle.orbit.velocity.vxy(idx), 'r', 'linewidth', 2);

% Draw the radial velocity
drawVector3d([particle.orbit.position.x(idx) particle.orbit.position.y(idx) 0 ], [particle.orbit.velocity.vrx(idx) particle.orbit.velocity.vry(idx) 0]/particle.orbit.velocity.vxy(idx), 'g', 'linewidth', 2);

xlabel('X axis')
ylabel('Y axis')
title(['Velocity XY at time ' num2str(particle.orbit.t(idx)*1E6) '\mu s'])
axis equal
view (2)
hold off
