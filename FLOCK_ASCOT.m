close all; clear all;
load test.mat
load('/home/andrea/ascot5/run_test.h5'); ascot = results.run_0517294857;
load('/home/andrea/Documents/MASTOrbit/AAAAAAAAAAAA_ORBIT_29976_at_-1.000_s_50_keV_60_deg_0_deg_R0.70_Z0.00.mat')

c = 299792458;
e = 1.602176462e-19;
%%%%%%%%%%%%% CALCULATING QUANTITIES FROM ASCOT  %%%%%%%%%%%%%
ascot.orbit.vnorm = sqrt(ascot.orbit.vr.^2 + ascot.orbit.vz.^2 + ascot.orbit.vphi.^2);
ascot.orbit.bnorm = sqrt(ascot.orbit.br.^2 + ascot.orbit.bz.^2 + ascot.orbit.bphi.^2);
ascot.orbit.pitch = (ascot.orbit.vr.*ascot.orbit.br + ascot.orbit.vz.*ascot.orbit.bz + ascot.orbit.vphi.*ascot.orbit.bphi)...
                     ./(ascot.orbit.vnorm.*ascot.orbit.bnorm);
ascot.orbit.x= ascot.orbit.r.*cos(ascot.orbit.phi.*pi./180);
ascot.orbit.y= ascot.orbit.r.*sin(ascot.orbit.phi.*pi./180);  
ascot.orbit.gamma = 1./ sqrt(1 -(ascot.orbit.vnorm.^2)./(c^2)); 
             
ascot.orbit.energy = (ascot.orbit.gamma-ones(1, length(ascot.orbit.gamma))).*...
(ascot.inistate.mass.*1.6605e-27).*(c.^2);
ascot.orbit.mu = ((ascot.inistate.mass.*1.6605e-27).* (ascot.orbit.gamma.*ascot.orbit.vnorm).^2)...
 .* (1 - ascot.orbit.pitch.^2)./(2.*ascot.orbit.bnorm);
ascot.orbit.vpar = ascot.orbit.pitch .*ascot.orbit.vnorm;




figure(1)
  subplot(3,3,1)
    plot(particle.orbit.t, particle.orbit.position.z, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.z, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('z position (m)')
    legend('FORCE', 'ASCOT')
    
  subplot(3,3,2)

    plot(particle.orbit.t, particle.orbit.position.R, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.r, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('r position (m)')
    legend('FORCE', 'ASCOT')
    

  subplot(3,3,3)
    plot(particle.orbit.t, particle.orbit.position.phi.*(180./pi), 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.phi, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('toroidal position (m)')
    legend('FORCE', 'ASCOT')
    
    subplot(3,3,4)
    plot(particle.orbit.t, particle.orbit.velocity.vz, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.vz, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('z speed (m/s)')
    legend('FORCE', 'ASCOT')
    
  subplot(3,3,5)

    plot(particle.orbit.t, particle.orbit.velocity.vr, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.vr, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('r speed (m/s)')
    legend('FORCE', 'ASCOT')
    

  subplot(3,3,6)
    plot(particle.orbit.t, particle.orbit.velocity.vtor, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.vphi, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('toroidal speed (m/s)')
    legend('FORCE', 'ASCOT')
    
    subplot(3,3,7)
    plot(particle.orbit.t, particle.orbit.position.x, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.x, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('x (m)')
    legend('FORCE', 'ASCOT')
    
  subplot(3,3,8)

   plot(particle.orbit.t, particle.orbit.position.y, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.y, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('y (m)')
    legend('FORCE', 'ASCOT')
    

  subplot(3,3,9)
    plot(particle.orbit.position.x, particle.orbit.position.y, 'b', 'linewidth', 2, ...
         ascot.orbit.x, ascot.orbit.y, 'r', 'linewidth', 2)
    xlabel('x (m)')
    ylabel('y (m)')
    legend('FORCE', 'ASCOT')

    
    
   figure(2)
   
 subplot(3,2,1)
    plot(particle.orbit.t, particle.orbit.magnetic_field.BR, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.br, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Br (Tesla)')
    legend('FORCE', 'ASCOT')
    
    subplot(3,2,2)
    plot(particle.orbit.t, particle.orbit.magnetic_field.BZ, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.bz, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Bz (Tesla)')
    legend('FORCE', 'ASCOT')
    
    subplot(3,2,3)
    plot(particle.orbit.t, particle.orbit.magnetic_field.BPHI, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.bphi, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Bphi (Tesla)')
    legend('FORCE', 'ASCOT')
    
     subplot(3,2,4)
    plot(particle.orbit.t,  particle.orbit.lambda, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.pitch, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('pitch')
    legend('FORCE', 'ASCOT')
    
    subplot(3,2,5)
    plot(particle.orbit.t,  particle.orbit.magnetic_moment.mu./particle.orbit.magnetic_moment.mu(1), 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.mu./ascot.orbit.mu(1), 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('\mu/\mu(1)')
    legend('FORCE', 'ASCOT')
    
    subplot(3,2,6)
    plot(particle.orbit.t, particle.orbit.energy.total./particle.orbit.energy.total(1), 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.energy./ascot.orbit.energy(1), 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('E_k/E_k(1)')
    legend('FORCE', 'ASCOT')
   
   figure(3)
   plot(particle.orbit.position.R, particle.orbit.position.z, 'b', 'linewidth', 2, ...
        ascot.orbit.r, ascot.orbit.z, 'r', 'linewidth', 2)
        xlabel('R (m)')
    ylabel('Z (m)')
        legend('FORCE', 'ASCOT')

   return
   

figure(3)
  subplot(2,3,1)
    plot(particle.orbit.t, particle.orbit.position.z, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.z, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('z position (m)')
    legend('FORCE', 'ASCOT')
    
  subplot(2,3,2)

    plot(particle.orbit.t, particle.orbit.position.R, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.r, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('r position (m)')
    legend('FORCE', 'ASCOT')
    

  subplot(2,3,3)
    plot(particle.orbit.t, particle.orbit.position.phi, 'b', 'linewidth', 2, ...
         ascot.orbit.time, ascot.orbit.phi./(180./pi), 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('toroidal position (m)')
    legend('FORCE', 'ASCOT')
    
    subplot(2,3,4)
    plot(particle.orbit.t, 100.*(1-particle.orbit.position.z./ascot.orbit.z), 'b', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('relative difference Z')
    legend('FORCE', 'ASCOT')
    
  subplot(2,3,5)

    plot(particle.orbit.t, 100.*(1-particle.orbit.position.R./ascot.orbit.r), 'b', 'linewidth', 2)

    xlabel('time (s)')
    ylabel('relative difference R')
    legend('FORCE', 'ASCOT')
    

  subplot(2,3,6)
    plot(particle.orbit.t, 100.*(1-particle.orbit.position.phi./(ascot.orbit.phi./(180./pi))), 'b', 'linewidth', 2)

    xlabel('time (s)')
    ylabel('relative difference phi')
    legend('FORCE', 'ASCOT')
   