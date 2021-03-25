function [particle] = MOrbit_gyro_average(particle)
% ------------------------------------------------------------------------
% GYRO-AVERARGE calculations
% ------------------------------------------------------------------------
% Determines how many gyro-orbits are done and their period and frequency
% particle.gyro_orbit.number = floor(t(end)*max(particle.orbit.angular_frequency));
% particle.gyro_orbit.number = abs(particle.gyro_orbit.theta/(2*pi));
tic();
printf('%s', 'Gyro-averaged calculations: '); fflush(stdout());
n = 1;
%phase = 2*pi*sign(diff(particle.gyro_orbit.theta_unwrapped(1:2)));
phase = 0;
for k = 1:particle.times.steps-1
  particle.gyro_orbit.theta(k) = particle.gyro_orbit.theta_unwrapped(k) - particle.gyro_orbit.theta_unwrapped(1);
  particle.gyro_orbit.ds(k) = sign(diff(particle.gyro_orbit.theta_unwrapped(k:k+1)));
  if (particle.gyro_orbit.ds(k) == 1 && particle.gyro_orbit.theta(k) > phase)
    particle.gyro_orbit.idx(n) = k;
    phase = phase + particle.gyro_orbit.ds(k)*2*pi;
    n = n + 1;
  endif
  if (particle.gyro_orbit.ds(k) == -1 && particle.gyro_orbit.theta(k) < phase)
    particle.gyro_orbit.idx(n) = k;
    phase = phase + particle.gyro_orbit.ds(k)*2*pi;
    n = n + 1;
  endif  
end  

if (n >= 2)

    particle.gyro_orbit.t = particle.orbit.t(particle.gyro_orbit.idx);
    particle.gyro_orbit.period = diff(particle.gyro_orbit.t);
    particle.gyro_orbit.number = length(particle.gyro_orbit.period);

    % Time array for the guiding centre calculations
    % particle.guidingcenter.t = particle.gyro_orbit.t(1:end-1) + particle.gyro_orbit.period/2;

    % Evaluates the guiding centre position and velocity
    for k = 1:particle.gyro_orbit.number
        particle.guidingcenter.t(k) = mean(particle.orbit.t(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));        
        particle.guidingcenter.position.x(k) = mean(particle.orbit.position.x(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.position.y(k) = mean(particle.orbit.position.y(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.position.z(k) = mean(particle.orbit.position.z(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.velocity.vx(k) = mean(particle.orbit.velocity.vx(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.velocity.vy(k) = mean(particle.orbit.velocity.vy(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.velocity.vz(k) = mean(particle.orbit.velocity.vz(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.velocity.v_par(k) = mean(particle.orbit.velocity.v_pll(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));
        particle.guidingcenter.velocity.v_per(k) = mean(particle.orbit.velocity.v_perp(particle.gyro_orbit.idx(k):particle.gyro_orbit.idx(k+1)));  
    end
    particle.guidingcenter.position.phi = atan2(particle.guidingcenter.position.y, particle.guidingcenter.position.x);
    particle.guidingcenter.position.phi = unwrap(particle.guidingcenter.position.phi);
    particle.guidingcenter.position.R = sqrt(particle.guidingcenter.position.x.^2 + particle.guidingcenter.position.y.^2);
    particle.guidingcenter.velocity.v = sqrt(particle.guidingcenter.velocity.vx.^2 +  particle.guidingcenter.velocity.vy.^2 + particle.guidingcenter.velocity.vz.^2);
    dt = particle.guidingcenter.t(2) - particle.guidingcenter.t(1);
    particle.guidingcenter.velocity.vtor = particle.guidingcenter.position.R.*gradient(particle.guidingcenter.position.phi, dt);

    % Evaluates the B field at the guiding centre position
    for k = 1:particle.gyro_orbit.number
        particle.guidingcenter.Bx(k) = Bfield(particle.guidingcenter.position.x(k), particle.guidingcenter.position.y(k), particle.guidingcenter.position.z(k))(1);
        particle.guidingcenter.By(k) = Bfield(particle.guidingcenter.position.x(k), particle.guidingcenter.position.y(k), particle.guidingcenter.position.z(k))(2);
        particle.guidingcenter.Bz(k) = Bfield(particle.guidingcenter.position.x(k), particle.guidingcenter.position.y(k), particle.guidingcenter.position.z(k))(3);  
    end
    particle.guidingcenter.B = sqrt(particle.guidingcenter.Bx.^2 + particle.guidingcenter.By.^2+particle.guidingcenter.Bz.^2);

    % Angular velocity
    particle.guidingcenter.angular_velocity = particle.charge*particle.guidingcenter.B/particle.mass;
    
    % Evaluates the averaged mu      
    particle.guidingcenter.larmor_radius = (particle.mass/abs(particle.charge))*particle.guidingcenter.velocity.v_per./particle.guidingcenter.B;
    particle.guidingcenter.magnetic_moment = particle.mass*(particle.guidingcenter.velocity.v_per.^2)./(2*particle.guidingcenter.B);
    particle.guidingcenter.gyroaveragetime = toc();
   
    
    % Compute pitch angle
    V1 = [particle.guidingcenter.velocity.vx', particle.guidingcenter.velocity.vy', particle.guidingcenter.velocity.vz'];
    V2 = [particle.guidingcenter.Bx', particle.guidingcenter.By', particle.guidingcenter.Bz'];
    particle.guidingcenter.pitch_angle = vectorAngle3d(V1, V2);
    clear V1 V2

    % Computes Lambda
    particle.guidingcenter.lambda = cos(particle.guidingcenter.pitch_angle);    

    printf('%f s\n', particle.guidingcenter.gyroaveragetime); fflush(stdout());

    % Evaluate the magnetic field gradient
    tic()
    printf('%s', 'Gyro-averaged grad B calculation: '); fflush(stdout());
    for k = 1:particle.gyro_orbit.number
    [gradB_x, gradB_y, gradB_z, gradBmod] = gradB(particle.guidingcenter.position.x(k), particle.guidingcenter.position.y(k), particle.guidingcenter.position.z(k));
    particle.guidingcenter.gradB_x(k) = gradB_x;
    particle.guidingcenter.gradB_y(k) = gradB_y;
    particle.guidingcenter.gradB_z(k) = gradB_z;
    particle.guidingcenter.gradBmod(k) = gradBmod;  
    end  

    clear gradB_x gradB_y gradB_z gradBmod

    % Calculates the relative variation of the magnetic field over a larmor radius
    particle.guidingcenter.dB = particle.guidingcenter.larmor_radius.*particle.guidingcenter.gradBmod;
    particle.guidingcenter.dB_over_B = particle.guidingcenter.dB./particle.guidingcenter.B;

    % Toroidal Canonical Angular Momentum
    particle.guidingcenter.PSI = PSI(particle.guidingcenter.position.R, particle.guidingcenter.position.z);
    particle.guidingcenter.Ptor = particle.mass*particle.guidingcenter.position.R.*particle.guidingcenter.velocity.vtor - particle.charge*particle.guidingcenter.PSI;
    
    particle.guidingcenter.gradCPUtime = toc();
    printf('%f s\n', particle.guidingcenter.gradCPUtime); fflush(stdout());

endif
