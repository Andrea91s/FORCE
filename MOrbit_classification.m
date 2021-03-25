function [particle] = MOrbit_classification(particle, EQLST, printout)
% Functions that classifies the orbits based on the following properties:
% 1) Lost or confined
% 2) Turning points (lambda)
% 3) Encircling


% Classifies a particle as confined or lost particles
[particle] = MOrbit_lost_particles(particle);
printf('... lost particles done!\n\n')

% It the particle is confined then check for the
% presence of turning points and encircling status
if (particle.orbit.classification.status1 == 0)
    % Classifies the particle as Bouncing or Not bouncing depending on
    % the presence of turning poitns (lambda = 0)
    [particle] = MOrbit_turning_points(particle);
    printf('... turning points done!\n\n')
endif
               

    % Classifies the particle as Encircling or Not Encircling depending
    % whether the orbit contains the magnetic axis or not
    [particle] = MOrbit_cyrcling_particles(particle, EQLST);
                printf('... cyrcling particles done!\n\n')
                
                

    % Classifies the particle in Passing, Stagnating, Trapped and Potato
    if      (particle.orbit.classification.status2 == 0 && particle.orbit.classification.status3 == 0)
            particle.orbit.classification.label = 'Stagnating';
            particle.orbit.classification.status = 2;   
    
    elseif  (particle.orbit.classification.status2 == 0 && particle.orbit.classification.status3 == 1)
            particle.orbit.classification.label = 'Passing';
            particle.orbit.classification.status = 1;
    
    elseif (particle.orbit.classification.status2 == 1 && particle.orbit.classification.status3 == 0)
            particle.orbit.classification.label = 'Banana';
            particle.orbit.classification.status = 3;   
    
    elseif (particle.orbit.classification.status2 == 1 && particle.orbit.classification.status3 == 1)
            particle.orbit.classification.label = 'Potato';
            particle.orbit.classification.status = 5;   
    
    endif


% Print the particle classification
if (printout == 1)
fid = 1;
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'ORBIT CLASSIFICATION\n')
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'Particle has status #1 = %d and is: %s\n', particle.orbit.classification.status1, particle.orbit.classification.label1);
fprintf(fid, 'Particle has status #2 = %d and is: %s\n', particle.orbit.classification.status2, particle.orbit.classification.label2);
fprintf(fid, 'Particle has status #3 = %d and is: %s\n', particle.orbit.classification.status3, particle.orbit.classification.label3);
fprintf(fid, '\n')
fprintf(fid, 'Particle has final status = %d and is: %s\n', particle.orbit.classification.status, particle.orbit.classification.label);
endif

