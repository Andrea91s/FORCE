directory = '/home/andrea/Documents/MASTORBIT/';
%orbit = 'ORBIT_29880_at_0.255_s_50_keV_37_deg_0_deg_R1.00_Z0.00.mat';
%orbit = 'ORBIT_29880_at_0.255_s_50_keV_66_deg_0_deg_R1.00_Z0.00.mat';
%orbit = 'ORBIT_29880_at_0.255_s_20_keV_37_deg_0_deg_R1.00_Z0.00.mat';
orbit = 'ORBIT_29880_at_0.255_s_50_keV_60_deg_R1.10_Z0.00.mat';

filename = [directory orbit];
load(filename);
MOrbit_fft(particle, particle.orbit.velocity.vx);

% Frequency analysis for the position
[F, P1x] = MOrbit_fft(particle, particle.orbit.position.x);
[F, P1y] = MOrbit_fft(particle, particle.orbit.position.y);
[F, P1z] = MOrbit_fft(particle, particle.orbit.position.z);

% Frequency analysis fot the velocity
[F, P1vx] = MOrbit_fft(particle, particle.orbit.velocity.vx);
[F, P1vy] = MOrbit_fft(particle, particle.orbit.velocity.vy);
[F, P1vz] = MOrbit_fft(particle, particle.orbit.velocity.vz);

% Plots
figure(1, 'position', [100 100 600 900]);
subplot(2,1,1)
loglog(F, [P1x P1y P1z]);
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')
legend('x', 'y', 'z')
title(strrep(orbit,'_', ' '))

subplot(2,1,2)
loglog(F, [P1vx P1vy P1vz]);
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')
legend('vx', 'vy', 'vz')

filename = ['ORBIT_' num2str(EQLST.pulseNumber) '_at_' num2str(EQLST.selectedtime, '%1.3f') '_s_' ...
                   num2str(particle.initialconditions.energy_eV/1000, '%2.0f') '_keV_' ...
                   num2str(180*particle.initialconditions.pitch_angle/pi, '%3.0f'), '_deg_' ...
                   'R' num2str(particle.orbit.position.R(1), '%1.2f') ...
                   '_Z' num2str(particle.orbit.position.z(1), '%1.2f') '.txt'];
       fid = fopen (filename, "w");