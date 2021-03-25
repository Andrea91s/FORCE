function [F, P1] = MOrbit_fft(particle, signal)
% function that evaluates the frequencies of the signal in input
% belonging to the particle object
%
% Example:
% directory = '/media/marco/WDpassport/marco/Documents/Orbit/MASTOrbit/Data/';
% filename = [directory 'ORBIT_29880_at_0.255_s_50_keV_37_deg_0_deg_R1.00_Z0.00.mat'];
% load(filename);
% MOrbit_fft(particle, particle.orbit.velocity.vx);


% Sampling frequency (in Hz)
Fs = 1/particle.times.dt;

% Signal to analyze
S = signal;

% Length of signal
L = length(S);

% FFT of the signal
Y = fft(S);

% Two-sided spectrum
P2 = abs(Y/L);

% Single-sided spectrum
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Frequencies in the one-sided spectrum (Hz)
F = Fs*(0:(L/2))/L;


% Makes the plot of the one-sided spectrum
loglog(F,P1)
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')
axis([min(F(2:end)) max(F) min(P1) max(P1)])


