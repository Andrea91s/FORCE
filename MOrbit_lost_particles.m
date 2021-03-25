function [particle] = MOrbit_lost_particles(particle)
% Function that check is the particle is lost by either hitting
% the coils and divertor or the wall

% Load the MAST geometry
MG = MOrbit_MAST_Geometry(fillyn = 1, 0);

% Check 
IN1U = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P1U(:,1), MG.Coils.P1U(:,2));
IN1L = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P1L(:,1), MG.Coils.P1L(:,2));

IN2U = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P2U(:,1), MG.Coils.P2U(:,2));
IN2L = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P2L(:,1), MG.Coils.P2L(:,2));

IN4U = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P4U(:,1), MG.Coils.P4U(:,2));
IN4L = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P4L(:,1), MG.Coils.P4L(:,2));

IN5U = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P5U(:,1), MG.Coils.P5U(:,2));
IN5L = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P5L(:,1), MG.Coils.P5L(:,2));

IN6U = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P6U(:,1), MG.Coils.P6U(:,2));
IN6L = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Coils.P6L(:,1), MG.Coils.P6L(:,2));

IND1 = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Divertor.part1(:,1), MG.Divertor.part1(:,2));
IND2 = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Divertor.part2(:,1), MG.Divertor.part2(:,2));
IND3 = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Divertor.part3(:,1), MG.Divertor.part3(:,2));
IND4 = inpolygon (particle.orbit.position.R, particle.orbit.position.z, MG.Divertor.part4(:,1), MG.Divertor.part4(:,2));

Vessel = [0.18 -2.2; 0.18 2.2; 2.0 2.2; 2.0 -2.2; 0.18 -2.2];
INVS = 1 - inpolygon (particle.orbit.position.R, particle.orbit.position.z, Vessel(:,1), Vessel(:,2));

IN = IN1U +  IN1L + IN2U + IN2L + IN4U + IN4L + IN5U + IN5L +IN6U + IN6L + IND1 + IND2 + IND3 + IND4 + INVS;

idx = find(IN > 0);
if (isempty(idx) == 1)
    particle.orbit.confinement = 'confined';
    particle.orbit.classification.label1 = 'Confined';
    particle.orbit.classification.status1 = 0;    
else
    particle.orbit.confinement = 'lost';
    particle.orbit.classification.label1 = 'Lost';
    particle.orbit.classification.status1 = 4;
    particle.orbit.classification.label2 = 'Lost';
    particle.orbit.classification.status2 = 4;
    particle.orbit.classification.label3 = 'Lost';
    particle.orbit.classification.status3 = 4;
    particle.orbit.classification.label = 'Lost';
    particle.orbit.classification.status = 4;
endif
