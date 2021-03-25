function [] = MOrbit_plot(particle, EQLST, two_screens)
% Function that plot the orbit of a particle

if (two_screens == 2)
  add_x = 1980;
else
  add_x = 0;
endif

% ------------------------------------------------------------------------------------------------------------------
% plot data vs time
% ------------------------------------------------------------------------------------------------------------------
figure(2, 'position', [100+add_x 100 1200 800])
  subplot(3,3,1)
    plot(particle.orbit.t, particle.orbit.position.x, 'b', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_x, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on
        plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.x, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.x, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif
    xlabel('time (s)')
    ylabel('X position (m)')
  
  subplot(3,3,4)
    plot(particle.orbit.t, particle.orbit.position.y, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_y, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.y, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.y, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif    
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif
    xlabel('time (s)')
    ylabel('Y position (m)')
    
  subplot(3,3,7)
    plot(particle.orbit.t, particle.orbit.position.z, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif    
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif
    xlabel('time (s)')
    ylabel('Z position (m)')
    

  subplot(3,3,2)
    plot(particle.orbit.t, particle.orbit.velocity.vx, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_vx, 'r', 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.vx, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif         
    xlabel('time (s)')
    ylabel('V_x (m/s)')
  
  subplot(3,3,5)
    plot(particle.orbit.t, particle.orbit.velocity.vy, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_vy, 'r', 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.vy, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif            
    xlabel('time (s)')
    ylabel('V_y (m/s)')
    
  subplot(3,3,8)
    plot(particle.orbit.t, particle.orbit.velocity.vz, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_vz, 'r', 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.vz, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif            
    xlabel('time (s)')
    ylabel('V_z (m/s)')    
    

  subplot(3,3,3)
    plot(particle.orbit.position.x, particle.orbit.position.y, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_y, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.x, particle.guidingcenter.position.y, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif       
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif  
    xlabel('X position (m)')
    ylabel('Y position (m)')
  
  subplot(3,3,6)
    plot(particle.orbit.position.x, particle.orbit.position.z, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.x, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif   
    xlabel('X position (m)')
    ylabel('Z position (m)')
    
  subplot(3,3,9)
    plot(particle.orbit.position.R, particle.orbit.position.z, 'linewidth', 2, ...
         particle.orbit.guidingcenter_R, particle.orbit.guidingcenter_z,  'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.R, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.R, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif   
    xlabel('Radial position (m)')
    ylabel('Z position (m)')    

    
    
    
% ------------------------------------------------------------------------------------------------------------------    
figure(3, 'position', [100+add_x 100 1200 800])
% ------------------------------------------------------------------------------------------------------------------
  subplot(3,3,1)
    plot(particle.orbit.t, particle.orbit.magnetic_field.Bx, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_Bx, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.Bx, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.Bx, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('B_x (Tesla)')
  
  subplot(3,3,4)
    plot(particle.orbit.t, particle.orbit.magnetic_field.By, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_By, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.By, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.By, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('B_y (Tesla)')
    
  subplot(3,3,7)
    plot(particle.orbit.t, particle.orbit.magnetic_field.Bz, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_Bz, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.Bz, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.Bz, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('B_z (Tesla)')
    

  subplot(3,3,2)
    plot(particle.orbit.position.x, particle.orbit.magnetic_field.B, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_B, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.B, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.x, particle.guidingcenter.B, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('X position (m)')
    ylabel('B (Tesla)')
  
  subplot(3,3,5)
    plot(particle.orbit.position.y, particle.orbit.magnetic_field.B, 'linewidth', 2, ...
        particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_B, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.y, particle.guidingcenter.gc.B, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.y, particle.guidingcenter.B, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('Y position (m)')
    ylabel('B (Tesla)')
    
  subplot(3,3,8)
    plot(particle.orbit.position.z, particle.orbit.magnetic_field.B, 'linewidth', 2, ...
         particle.orbit.guidingcenter_z, particle.orbit.guidingcenter_B, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.z, particle.guidingcenter.gc.B, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.z, particle.guidingcenter.B, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif     
    xlabel('z position (m)')
    ylabel('B (Tesla)')

  subplot(3,3,3)
    plot(particle.orbit.t, particle.orbit.magnetic_field.B, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_B, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.B, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.B, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif     
    xlabel('time (s)')
    ylabel('B (Tesla)')
  
  subplot(3,3,6)
    plot(particle.orbit.t, particle.orbit.larmor_radius, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_rl, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(    
	 particle.guidingcenter.gc.t, particle.guidingcenter.gc.larmor_radius, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.larmor_radius, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif     
    xlabel('time (s)')
    ylabel('Larmor Radius (m)')
    
  subplot(3,3,9)
    plot(particle.orbit.t, particle.orbit.angular_velocity, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_omega, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(
	 particle.guidingcenter.gc.t, particle.guidingcenter.gc.angular_velocity, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.angular_velocity, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif     
    xlabel('time (s)')
    ylabel('\omega (rad/s)')  
    
    
    
% ------------------------------------------------------------------------------------------------------------------
figure(4, 'position', [100+add_x 100 1200 800])
% ------------------------------------------------------------------------------------------------------------------
  subplot(3,3,1)
    plot(particle.orbit.t, particle.orbit.magnetic_field.B, 'linewidth', 2, 
         particle.orbit.t, particle.orbit.guidingcenter_B, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.B, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.B, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('B (Tesla)')
  
  subplot(3,3,4)
    plot(particle.orbit.t, particle.orbit.magnetic_moment.mu/physical_constant ('electron volt'), 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_mu/physical_constant ('electron volt'), 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on             
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.magnetic_moment.mu/physical_constant ('electron volt'), 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.magnetic_moment/physical_constant ('electron volt'), 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('\mu (eV/Tesla)')  
    
  subplot(3,3,7)
    plot(particle.orbit.t, particle.orbit.larmor_radius, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_rl, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.larmor_radius, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.larmor_radius, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('Larmor Radius (m)')
    

  subplot(3,3,2)
    plot(particle.orbit.t, particle.orbit.velocity.v_pll, 'b', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.velocity.v_perp, 'b--', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter.v_parl, 'r', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter.v_perp, 'r--', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.velocity.v_per, 'g--', 'linewidth', 2, ...
              particle.guidingcenter.gc.t, particle.guidingcenter.gc.velocity.v_par, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.v_par, 'color', [0. 0.5 0.], 'linewidth', 2, ...
                 particle.guidingcenter.t, particle.guidingcenter.velocity.v_per, 'color', [0. 0.5 0.], 'linewidth', 2, 'linestyle', '--')
        hold off
    endif      
    xlabel('time (s)')
    ylabel('V (m/s)')
    legend('v_{pll}', 'v_{perp}')
  
  subplot(3,3,5)
    plot(particle.orbit.t, particle.orbit.lambda, 'linewidth', 2, ...
	 particle.orbit.t, particle.orbit.guidingcenter.lambda, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.lambda, 'g', 'linewidth', 2)
        hold off
    endif
%    if (isfield(particle.gyro_orbit, 't') == 1)
%        hold on
%            plot(particle.guidingcenter.t, particle.guidingcenter.v_par, 'color', [0. 0.5 0.], 'linewidth', 2, ...
%                 particle.guidingcenter.t, particle.guidingcenter.v_per, 'color', [0. 0.5 0.], 'linewidth', 2, 'linestyle', '--')
%        hold off
%    endif      
    xlabel('time (s)')
    ylabel('Lambda')
    
  subplot(3,3,8)
    plot(particle.orbit.t, particle.orbit.pitch_angle, 'linewidth', 2, ...
    	 particle.orbit.t, particle.orbit.guidingcenter.pitch_angle, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.pitch_angle, 'g', 'linewidth', 2)
        hold off
    endif
    xlabel('time (s)')
    ylabel('Pitch angle (rad)')  

  subplot(3,3,3)
    plot(particle.orbit.t, particle.orbit.energy.kinetic/physical_constant ('electron volt'), 'r','linewidth', 2, ...
         particle.orbit.t, particle.orbit.energy.potential/physical_constant ('electron volt'), 'b', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.energy.total/physical_constant ('electron volt'), 'k', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Total Energy (eV)')
    legend('E_{Kin}', 'E_{Pot}', 'E_{Tot}')
  
  subplot(3,3,6)
    plot(particle.orbit.t, particle.orbit.energy.kinetic_perp/physical_constant ('electron volt'), 'b--','linewidth', 2, ...
         particle.orbit.t, particle.orbit.energy.kinetic_parl/physical_constant ('electron volt'), 'b', 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.energy.kinetic/physical_constant ('electron volt'), 'm', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.energy.kinetic_per/physical_constant ('electron volt'), 'g--', 'linewidth', 2, ...
         particle.guidingcenter.gc.t, particle.guidingcenter.gc.energy.kinetic_par/physical_constant ('electron volt'), 'g', 'linewidth', 2, ...
         particle.guidingcenter.gc.t, particle.guidingcenter.gc.energy.kinetic/physical_constant ('electron volt'), 'k-.', 'linewidth', 2)
        hold off
    endif
    xlabel('time (s)')
    ylabel('Kinetic Energy (eV)')
    legend('E_{PERP}', 'E_{PARL}', 'E_{Tot}')
    
  subplot(3,3,9)
    %plot(particle.orbit.t, particle.orbit.magnetic_field.dB_over_B, 'linewidth', 2)
    %plot(particle.guidingcenter.t, particle.guidingcenter.dB_over_B, 'linewidth', 2, ...
     plot(particle.orbit.t, particle.orbit.guidingcenter_deltaBoverB, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('\delta B / B')  
    
    
    
    
% ------------------------------------------------------------------------------------------------------------------ 
% Plot the particle trajectory on the plane perpendicular to B
% moving with the guiding centre at its origin
% ------------------------------------------------------------------------------------------------------------------
figure(5, 'position', [100+add_x 100 600 480])
  plot(particle.orbit.larmor_radius.*cos(particle.gyro_orbit.theta_unwrapped ), particle.orbit.larmor_radius.*sin(particle.gyro_orbit.theta_unwrapped), 'linewidth', 2, ...
       particle.orbit.larmor_radius(1)*cos(particle.gyro_orbit.theta_unwrapped(1)), particle.orbit.larmor_radius(1)*sin(particle.gyro_orbit.theta_unwrapped(1)), 'ob', 'linewidth', 2, ...
       particle.orbit.larmor_radius(end)*cos(particle.gyro_orbit.theta_unwrapped(end)), particle.orbit.larmor_radius(end)*sin(particle.gyro_orbit.theta_unwrapped(end)), 'or', 'linewidth', 2)   
title('Projection of the orbit on the plane perpendicular to B')    
    
    

% ------------------------------------------------------------------------------------------------------------------
% Plot the poloidal projection of the trajectory
% ------------------------------------------------------------------------------------------------------------------
figure(6, 'position', [100+add_x 100 600 480])
 plot(particle.orbit.position.R, particle.orbit.position.z, 'b', 'linewidth', 2, ...
      particle.orbit.guidingcenter_R, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2, ...
      EQLST.magneticAxisr, EQLST.magneticAxisz, 'k+', 'markersize', 12)
    if (particle.guidingcenter.gc.solve == 1)    
	hold on         
         %plot(particle.guidingcenter.gc.position.R, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    %if (isfield(particle.gyro_orbit, 't') == 1)
     %   hold on
      %      plot(particle.guidingcenter.position.R, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
       % hold off
    %endif      
%if (isfield(isfield(particle.orbit, 'turning_points') == 1 && particle.orbit.turning_points, 'x') == 1)
%       hold on
%	    plot(particle.orbit.turning_points.t, particle.orbit.turning_points.x, 'k+', 'markersize', 10, 'linewidth', 2)
%	 hold off
%    endif     
%axis equal
xlabel('R (m)')
ylabel('Z (m)')
MG = MAST_Geometry(fillyn = 1);
hold on
 contour(EQLST.R, EQLST.Z, EQLST.poloidalFlux, 20, 'k');
 plot(EQLST.rb, EQLST.zb, 'k--', 'linewidth', 2) 
hold off
title(['Poloidal Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])




% ------------------------------------------------------------------------------------------------------------------
% Plot the equatorial projection of the trajectory
% ------------------------------------------------------------------------------------------------------------------
figure(7,  'position', [100+add_x 100 600 480])
    plot(particle.orbit.position.x, particle.orbit.position.y, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_y, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.position.x, particle.guidingcenter.position.y, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('X position (m)')
    ylabel('Y position (m)')
    hold on
      MAST_Geometry_Equatorial()
    hold off
    title(['Equatorial Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])




   
% ------------------------------------------------------------------------------------------------------------------    
% Plot the normalized variation of the constant of motion    
% ------------------------------------------------------------------------------------------------------------------
figure(8, 'position', [100+add_x 100 900 400])

  subplot(1,2,1)
    plot(particle.orbit.t, particle.orbit.magnetic_moment.mu/particle.orbit.magnetic_moment.mu(1), 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_mu/particle.orbit.guidingcenter_mu(1), 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.magnetic_moment.mu/particle.guidingcenter.gc.magnetic_moment.mu(1), 'g', 'linewidth', 2)
        hold off
    endif
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.magnetic_moment/particle.guidingcenter.magnetic_moment(1), 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('\mu/\mu(0) (eV/Tesla)')  
    

  subplot(1,2,2)
    plot(particle.orbit.t, particle.orbit.energy.kinetic/particle.orbit.energy.kinetic(1), 'b', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.energy.kinetic/particle.guidingcenter.gc.energy.kinetic(1), 'r', 'linewidth', 2)
        hold off
    endif         
    xlabel('time (s)')
    ylabel('Total Kinetic Energy /E_{K,tot}(0) (eV)')
    

    
   
% ------------------------------------------------------------------------------------------------------------------    
% Plot the toroidal velocity and angle   
% ------------------------------------------------------------------------------------------------------------------
figure(9, 'position', [100+add_x 100 900 400])

  subplot(1,2,1)
    plot(particle.orbit.t, particle.orbit.velocity.vtor, 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.vtor, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('Toroidal Velocity (m/s)')  
    

  subplot(1,2,2)
    plot(particle.orbit.t, particle.orbit.position.phi, 'b', 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.phi, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif          
    xlabel('time (s)')
    ylabel('Toroidal Angle (rad)')  


% ------------------------------------------------------------------------------------------------------------------
% Guiding centre approximation plots
% ------------------------------------------------------------------------------------------------------------------
if (particle.guidingcenter.gc.solve == 1)
% plot data vs time
figure(10, 'position', [100+add_x 100 1200 800])
  subplot(3,3,1)
    plot(particle.orbit.t, particle.errror_x, 'k', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Error in X position (m)')
  
  subplot(3,3,4)
    plot(particle.orbit.t, particle.errror_y, 'k', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Error in Y position (m)')
    
  subplot(3,3,7)
    plot(particle.orbit.t, particle.errror_z, 'k', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Error in X position (m)')
    xlabel('time (s)')
    ylabel('Error in Z position (m)')
    

  subplot(3,3,2)
    plot(particle.orbit.t, particle.orbit.guidingcenter_V, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.velocity.vmod, 'g', 'linewidth', 2)
        hold off
    endif
    xlabel('time (s)')
    ylabel('GC velocity (m/s)')
  
  subplot(3,3,5)
    plot(particle.orbit.t, particle.orbit.velocity.vy, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_vy, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('V_y (m/s)')
    
  subplot(3,3,8)
    plot(particle.orbit.t, particle.orbit.velocity.vz, 'linewidth', 2, ...
         particle.orbit.t, particle.orbit.guidingcenter_vz, 'r', 'linewidth', 2)
    xlabel('time (s)')
    ylabel('V_z (m/s)')    
    

  subplot(3,3,3)
    plot(particle.orbit.position.x, particle.orbit.position.y, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_y, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, 'g', 'linewidth', 2)
        hold off
    endif
    xlabel('X position (m)')
    ylabel('Y position (m)')
  
  subplot(3,3,6)
    plot(particle.orbit.position.x, particle.orbit.position.z, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    xlabel('X position (m)')
    ylabel('Z position (m)')
    
  subplot(3,3,9)
    plot(particle.orbit.position.y, particle.orbit.position.z, 'linewidth', 2, ...
         particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_z,  'r', 'linewidth', 2)
    if (particle.guidingcenter.gc.solve == 1)
	hold on         
         plot(
         particle.guidingcenter.gc.position.y, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
        hold off
    endif
    xlabel('Y position (m)')
    ylabel('Z position (m)')    


endif




% ------------------------------------------------------------------------------------------------------------------
% Plot the 3D trajectory within MAST vessel
% ------------------------------------------------------------------------------------------------------------------
figure(20, 'position', [100+add_x 100 600 480])
plot3(particle.orbit.position.x, particle.orbit.position.y, particle.orbit.position.z, 'b', 'linewidth', 3)
hold on
  plot3(particle.orbit.position.x(1), particle.orbit.position.y(1), particle.orbit.position.z(1), 'ob', 'linewidth', 2)
  plot3(particle.orbit.position.x(end), particle.orbit.position.y(end), particle.orbit.position.z(end), 'or', 'linewidth', 2)
  plot3(particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_y, particle.orbit.guidingcenter_z, 'r', 'linewidth',3);
  if (particle.guidingcenter.gc.solve == 1)
    plot3(particle.guidingcenter.gc.position.x, particle.guidingcenter.gc.position.y, particle.guidingcenter.gc.position.z, 'g', 'linewidth', 2)
  endif
    if (isfield(particle.gyro_orbit, 't') == 1)
            plot3(particle.guidingcenter.position.x, particle.guidingcenter.position.y, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
    endif      
%  if (isfield(particle.orbit, 'turning_points') == 1 && isfield(particle.orbit.turning_points, 'x') == 1)
%    plot3(particle.orbit.turning_points.x, particle.orbit.turning_points.y, particle.orbit.turning_points.z, 'bo',  'markerfacecolor', 'b', 'markersize', 10, 'linewidth', 2)
%    %plot3(particle.orbit.turning_points.gc_x, particle.orbit.turning_points.gc_y, particle.orbit.turning_points.gc_z, 'ro',  'markerfacecolor', 'r', 'markersize', 10, 'linewidth', 2)
%  endif
%  if (particle.guidingcenter.gc.solve == 1)
%    if (isfield(particle.orbit, 'turning_points') == 1 && isfield(particle.orbit.turning_points, 'gc') == 1)
%      plot3(particle.orbit.turning_points.gc.x, particle.orbit.turning_points.gc.y, particle.orbit.turning_points.gc.z, 'go',  'markerfacecolor', 'g', 'markersize', 10, 'linewidth', 2)
%    endif
%  endif
  plot3(particle.orbit.guidingcenter_R, 0*particle.orbit.guidingcenter_R, particle.orbit.guidingcenter_z, 'g', 'linewidth', 3)
    if (isfield(particle.gyro_orbit, 't') == 1)
            plot3(particle.guidingcenter.position.R, 0*particle.guidingcenter.position.R, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
    endif    
  %drawCylinder([0 0 -1.7 0 0 1.7 0.15],'facecolor', [0.9 0.9 0.9])

  
  % Draws the last closed flux surface
  [x y z] = revolutionSurface([EQLST.rb' EQLST.zb'], linspace(0, pi, 90));
  LCFScolor = [227 174 204]/256;
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  %set(h, 'facecolor', 'none')
  set(h, 'facecolor', LCFScolor)
  
  % Draws the polidal field coils
  coilcolor = [0.914 0.902 0.278];
  [x y z] = revolutionSurface(MG.Coils.P2U, linspace(0, 2*pi, 180));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P4U, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P5U, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)  
    [x y z] = revolutionSurface(MG.Coils.P6U, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P2L, linspace(0, 2*pi, 180));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P4L, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
  [x y z] = revolutionSurface(MG.Coils.P5L, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)  
    [x y z] = revolutionSurface(MG.Coils.P6L, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', coilcolor)
      [x y z] = revolutionSurface(  MG.Divertor.part1, linspace(0, pi, 90));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', [0.8 0.8 0.8])  
      [x y z] = revolutionSurface(  MG.Divertor.part2, linspace(0, 2*pi, 180));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', [0.8 0.8 0.8])    
      [x y z] = revolutionSurface(  MG.Divertor.part3, linspace(0, 2*pi, 180));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', [0.8 0.8 0.8])  
      [x y z] = revolutionSurface(  MG.Divertor.part4, linspace(0, 2*pi, 180));
  h = surf(x,y,z);
  set(h, 'linewidth', 0.5, 'edgecolor', 'k')
  set(h, 'facecolor', [0.8 0.8 0.8])      
hold off
axis equal
axis off
box off
hidden
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Z position (m)')
view(308.06, 12.883)

title(['Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])


% ------------------------------------------------------------------------------------------------------------------
% Plot the 3D normalized vectors of the B field and of the velocity
% ------------------------------------------------------------------------------------------------------------------
figure(50); 
hold on
hB = drawVector3d ([0 0 0], [particle.initialconditions.bx particle.initialconditions.by particle.initialconditions.bz]);
set(hB, 'linewidth', 2)
hV = drawVector3d ([0 0 0], [particle.initialconditions.vx0/particle.initialconditions.v particle.initialconditions.vy0/particle.initialconditions.v particle.initialconditions.vz0/particle.initialconditions.v]);
set(hV, 'linewidth', 2, 'color', 'r')
xlabel('x')
ylabel('y')
zlabel('z')
title('Initial pitch-angle')
legend('B field', 'Velocity')
hold off




endfunction
