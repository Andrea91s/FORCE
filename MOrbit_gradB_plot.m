function [] = MOrbit_gradB_plot(particle)
figure(60); 
hold on
hB = drawVector3d ([0 0 0], [particle.initialconditions.bx particle.initialconditions.by particle.initialconditions.bz]);
set(hB, 'linewidth', 2)
hV = drawVector3d ([0 0 0], [particle.initialconditions.vx0/particle.initialconditions.v particle.initialconditions.vy0/particle.initialconditions.v particle.initialconditions.vz0/particle.initialconditions.v]);
set(hV, 'linewidth', 2, 'color', 'r')
hG = drawVector3d ([0 0 0], [particle.orbit.magnetic_field.gradB_x(1)/particle.orbit.magnetic_field.gradBmod(1) ...
                             particle.orbit.magnetic_field.gradB_y(1)/particle.orbit.magnetic_field.gradBmod(1) ...
                             particle.orbit.magnetic_field.gradB_z(1)/particle.orbit.magnetic_field.gradBmod(1)]);
set(hG, 'linewidth', 2, 'color', 'g')
xlabel('x')
ylabel('y')
zlabel('z')
title('Initial pitch-angle')
legend('B field', 'Velocity', '\nabla B')
hold off

