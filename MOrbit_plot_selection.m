function [] = MOrbit_plot_selection(particle, EQLST, two_screens)
% Function that plot only those data necessary to determine if the
% the orbit is passing or trapped

if (two_screens == 2)
  add_x = 1980;
else
  add_x = 0;
endif

 figure(1)
subplot(1, 2, 1) 
 plot(particle.orbit.position.R, particle.orbit.position.z, 'b', 'linewidth', 2, ...
         particle.orbit.guidingcenter_R, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2, ...
         EQLST.magneticAxisr, EQLST.magneticAxisz, 'm+', 'markersize', 12, 'linewidth', 2)
    %if (isfield(particle.gyro_orbit, 't') == 1)
    %    hold on
    %        plot(particle.guidingcenter.position.R, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
    %    hold off
    %endif         
    xlabel('R (m)')
    ylabel('Z (m)')
    MG = MAST_Geometry(fillyn = 1);
    hold on
    contour(EQLST.R, EQLST.Z, EQLST.poloidalFlux, 20, 'k');
    plot(EQLST.rb, EQLST.zb, 'k--', 'linewidth', 2) 
    hold off
    title(['Poloidal Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])

subplot(1, 2, 2)  
plot(particle.orbit.t, particle.orbit.position.phi, 'b', 'linewidth', 2)
        if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.phi, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif 
    xlabel('time (s)')
    ylabel('toroidal angle \phi (rad)')
    title(['Toroidal Angle - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.']) 
return


 
    
% ------------------------------------------------------------------------------------------------------------------    
figure(30, 'position', [50+add_x 100 1800 800])
% ------------------------------------------------------------------------------------------------------------------


  subplot(2,4,1)
    plot(particle.orbit.t, particle.orbit.lambda, 'linewidth', 2)
    xlabel('time (s)')
    ylabel('Lambda')



  subplot(2,4,2)
    plot(particle.orbit.t, particle.orbit.magnetic_moment.mu/particle.orbit.magnetic_moment.mu(1), 'linewidth', 2)
    %if (isfield(particle.gyro_orbit, 't') == 1)
        %hold on
        %    plot(particle.guidingcenter.t, particle.guidingcenter.magnetic_moment/particle.guidingcenter.magnetic_moment(1), 'color', [0. 0.5 0.], 'linewidth', 2)
        %hold off
    %endif  
    axis([particle.orbit.t(1) particle.orbit.t(end)])
    h = xlabel('time (s)'); set(h, 'fontsize', 12)
    h = ylabel('\mu (t)/\mu(0)'); set(h, 'fontsize', 12)  
    
    

  subplot(2,4,3)
    plot(particle.orbit.t, particle.orbit.energy.kinetic/particle.orbit.energy.kinetic(1), 'b', 'linewidth', 2)      
        axis([particle.orbit.t(1) particle.orbit.t(end)])
    h = xlabel('time (s)'); set(h, 'fontsize', 12)
    h = ylabel('Total Kinetic Energy E_{K}(t)/E_{K,tot}(0)'); set(h, 'fontsize', 12)
 

  subplot(2,4,4)
    plot(particle.orbit.t, particle.orbit.Ptor/particle.orbit.Ptor(1), 'b', 'linewidth', 2)  
        axis([particle.orbit.t(1) particle.orbit.t(end)])
    h = xlabel('time (s)'); set(h, 'fontsize', 12)
    h = ylabel('P_{\phi} (t) /P_{\phi}(0)') ; set(h, 'fontsize', 12)

    
    
  subplot(2,4,5)
    plot(particle.orbit.position.R, particle.orbit.position.z, 'b', 'linewidth', 2, ...
         particle.orbit.guidingcenter_R, particle.orbit.guidingcenter_z, 'r', 'linewidth', 2, ...
         EQLST.magneticAxisr, EQLST.magneticAxisz, 'm+', 'markersize', 12, 'linewidth', 2)
    %if (isfield(particle.gyro_orbit, 't') == 1)
    %    hold on
    %        plot(particle.guidingcenter.position.R, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
    %    hold off
    %endif         
    xlabel('R (m)')
    ylabel('Z (m)')
    MG = MAST_Geometry(fillyn = 1);
    hold on
    contour(EQLST.R, EQLST.Z, EQLST.poloidalFlux, 20, 'k');
    plot(EQLST.rb, EQLST.zb, 'k--', 'linewidth', 2) 
    hold off
    title(['Poloidal Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])




  subplot(2,4,6)
    plot(particle.orbit.position.x, particle.orbit.position.y, 'linewidth', 2, ...
         particle.orbit.guidingcenter_x, particle.orbit.guidingcenter_y, 'r', 'linewidth', 2)
    %if (isfield(particle.gyro_orbit, 't') == 1)
    %hold on
    %        plot3(particle.guidingcenter.position.x, particle.guidingcenter.position.y, particle.guidingcenter.position.z, 'color', [0. 0.5 0.], 'linewidth', 2)
    %hold off
    %endif   
    xlabel('X position (m)')
    ylabel('Y position (m)')
    hold on
      MAST_Geometry_Equatorial()
    hold off
    title(['Equatorial Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])

    
  subplot(2,4,7)
    plot(particle.orbit.t, particle.orbit.position.phi, 'b', 'linewidth', 2)
        if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.position.phi, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif 
    xlabel('time (s)')
    ylabel('toroidal angle \phi (rad)')
    title(['Toroidal Angle - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])    

  subplot(2,4,8)
    plot(particle.orbit.t, particle.orbit.velocity.vtor, 'linewidth', 2)
    if (isfield(particle.gyro_orbit, 't') == 1)
        hold on
            plot(particle.guidingcenter.t, particle.guidingcenter.velocity.vtor, 'color', [0. 0.5 0.], 'linewidth', 2)
        hold off
    endif      
    xlabel('time (s)')
    ylabel('Toroidal Velocity (m/s)')  
    


%  subplot(2,3,2)
%    plot(particle.orbit.t, particle.orbit.pitch_angle, 'linewidth', 2, ...
%    	 particle.orbit.t, particle.orbit.guidingcenter.pitch_angle, 'r', 'linewidth', 2)
%    if (particle.guidingcenter.gc.solve == 1)
%	hold on         
%         plot(particle.guidingcenter.gc.t, particle.guidingcenter.gc.pitch_angle, 'g', 'linewidth', 2)
%        hold off
%    endif
%    xlabel('time (s)')
%    ylabel('Pitch angle (rad)')  
