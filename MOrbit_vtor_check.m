function [particle] = MOrbit_vtor_check(particle)
% Function that tries to determine if the toroidal velocity of a particle
% changes sign indicating the presence of turning points

tic();
printf('%s', 'Toroidal velocity turning points: '); fflush(stdout());


% Check the condition for the instantaneuous toroidal velocity
particle.orbit.vtor_reflection.t = zerocrossing (particle.orbit.t, particle.orbit.velocity.vtor);

if (isempty(particle.orbit.vtor_reflection.t) != 1)
  particle.orbit.vtor_reflection.status = 'Instantaneuous reflection';
  particle.orbit.vtor_reflection.N = length(particle.orbit.vtor_reflection.t);
  particle.orbit.vtor_reflection.bounces = particle.orbit.vtor_reflection.N  - 1;


  % Find the indexes in the time array where these are:
  for n = 1:particle.orbit.vtor_reflection.N
    particle.orbit.vtor_reflection.idx(n) = max(find(particle.orbit.t <= particle.orbit.vtor_reflection.t(n)));
  end  
  particle.orbit.vtor_reflection.x = interp1(particle.orbit.t, particle.orbit.position.x, particle.orbit.vtor_reflection.t);
  particle.orbit.vtor_reflection.y = interp1(particle.orbit.t, particle.orbit.position.y, particle.orbit.vtor_reflection.t);
  particle.orbit.vtor_reflection.z = interp1(particle.orbit.t, particle.orbit.position.z, particle.orbit.vtor_reflection.t);
  particle.orbit.vtor_reflection.R = sqrt(particle.orbit.vtor_reflection.x.^2 + particle.orbit.vtor_reflection.y.^2);

  particle.orbit.vtor_reflection.bounce_Bmax = interp1(particle.orbit.t, particle.orbit.magnetic_field.B, particle.orbit.vtor_reflection.t);
  particle.orbit.vtor_reflection.bounce_Bmin(1) = min(particle.orbit.magnetic_field.B(1:particle.orbit.vtor_reflection.idx(1)));
  if (particle.orbit.vtor_reflection.N > 1)
     for n = 2:particle.orbit.vtor_reflection.N
         particle.orbit.vtor_reflection.bounce_Bmin(n) = min(particle.orbit.magnetic_field.B(particle.orbit.vtor_reflection.idx(n-1):particle.orbit.vtor_reflection.idx(n)));
     end
  end
  else
    particle.orbit.vtor_reflection.status = 'No instantaneuous reflection';
endif





% Check the condition for the gyro-averaged toroidal velocity
particle.guidingcenter.vtor_reflection.t = zerocrossing (particle.guidingcenter.t, particle.guidingcenter.velocity.vtor);

if (isempty(particle.guidingcenter.vtor_reflection.t) != 1)
  particle.guidingcenter.vtor_reflection.status = 'Gyro-average reflection';
  particle.guidingcenter.vtor_reflection.N = length(particle.guidingcenter.vtor_reflection.t);
  particle.guidingcenter.vtor_reflection.bounces = particle.guidingcenter.vtor_reflection.N  - 1;


  % Find the indexes in the time array where these are:
  for n = 1:particle.guidingcenter.vtor_reflection.N
    particle.guidingcenter.vtor_reflection.idx(n) = max(find(particle.orbit.t <= particle.guidingcenter.vtor_reflection.t(n)));
  end  
  particle.guidingcenter.vtor_reflection.x = interp1(particle.guidingcenter.t, particle.guidingcenter.position.x, particle.guidingcenter.vtor_reflection.t);
  particle.guidingcenter.vtor_reflection.y = interp1(particle.guidingcenter.t, particle.guidingcenter.position.y, particle.guidingcenter.vtor_reflection.t);
  particle.guidingcenter.vtor_reflection.z = interp1(particle.guidingcenter.t, particle.guidingcenter.position.z, particle.guidingcenter.vtor_reflection.t);
  particle.guidingcenter.vtor_reflection.R = sqrt(particle.guidingcenter.vtor_reflection.x.^2 + particle.guidingcenter.vtor_reflection.y.^2);

  particle.guidingcenter.vtor_reflection.bounce_Bmax = interp1(particle.orbit.t, particle.orbit.magnetic_field.B, particle.guidingcenter.vtor_reflection.t);
  particle.guidingcenter.vtor_reflection.bounce_Bmin(1) = min(particle.orbit.magnetic_field.B(1:particle.guidingcenter.vtor_reflection.idx(1)));
  if (particle.guidingcenter.vtor_reflection.N > 1)
     for n = 2:particle.guidingcenter.vtor_reflection.N
         particle.guidingcenter.vtor_reflection.bounce_Bmin(n) = min(particle.orbit.magnetic_field.B(particle.guidingcenter.vtor_reflection.idx(n-1):particle.guidingcenter.vtor_reflection.idx(n)));
     end
  end
  else
    particle.guidingcenter.vtor_reflection.status = 'No gyro-average reflection';
endif

printf('%f s\n', toc()); fflush(stdout());
