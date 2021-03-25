function [] = MOrbit_plot_selection(particle, EQLST, two_screens)
% Function that plot only those data necessary to determine if the
% the orbit is passing or trapped

if (two_screens == 2)
  add_x = 1980;
else
  add_x = 0;
endif

   close all 
    
    
% ------------------------------------------------------------------------------------------------------------------    
figure(30, 'position', [50+add_x 100 1800 800])
% ------------------------------------------------------------------------------------------------------------------

    
       plot(EQLST.magneticAxisr, EQLST.magneticAxisz, 'k+', 'markersize', 25, 'linewidth', 3)
    xlabel('R (m)')
    ylabel('Z (m)')
    MG = MAST_Geometry(fillyn = 1);
    hold on
    contour(EQLST.R, EQLST.Z, EQLST.poloidalFlux, 20, 'k');
    plot(EQLST.rb, EQLST.zb, 'r', 'linewidth', 3) 
    hold off
    title(['Poloidal Projection - Pulse ' num2str(EQLST.pulseNumber)  ', at ' num2str(EQLST.selectedtime), ' s.'])



    
  