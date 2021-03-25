function [particle] = MOrbit_gradB(particle)
% Calculate the gradients along the particle trajectory

    tic()
    printf('%s', 'Grad B calcultion: ');
    for t = 1:particle.times.steps
        [gradB_x, gradB_y, gradB_z, gradBmod] = gradB(particle.orbit.position.x(t),particle.orbit.position.y(t),particle.orbit.position.z(t));
        particle.orbit.magnetic_field.gradB_x(t) = gradB_x;
        particle.orbit.magnetic_field.gradB_y(t) = gradB_y;  
        particle.orbit.magnetic_field.gradB_z(t) = gradB_z;
        particle.orbit.magnetic_field.gradBmod(t) = gradBmod;
    end 
    particle.orbit.gradCPUtime = toc();
    clear gradB_x gradB_y gradB_z gradBmod

    % Calculates the relative variation of the magnetic field over a larmor radius
    particle.orbit.magnetic_field.dB = particle.orbit.larmor_radius.*particle.orbit.magnetic_field.gradBmod';
    particle.orbit.magnetic_field.dB_over_B = particle.orbit.magnetic_field.dB./particle.orbit.magnetic_field.B;    
    printf('%f s\n', toc());
