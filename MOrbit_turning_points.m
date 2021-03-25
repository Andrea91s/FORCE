function [particle] = MOrbit_turning_points(particle)
% ------------------------------------------------------------------------
% Turning points location, bounce frequency and mirror ratio
% ------------------------------------------------------------------------
tic();
printf('%s', 'Turning points: '); fflush(stdout());
% ORBIT
% Search for where lambda = 0, that is when the parellel
% velocity is zero: lambda = cos(vpll/v)
particle.orbit.turning_points.t = zerocrossing (particle.orbit.t, particle.orbit.lambda);

% Calculate the location of the turning points
if (isempty(particle.orbit.turning_points.t) != 1)

  particle.orbit.turning_points.status = 'trapped';
  particle.orbit.classification.label2 = 'Bouncing';
  particle.orbit.classification.status2 = 1;
  particle.orbit.turning_points.N = length(particle.orbit.turning_points.t);
  particle.orbit.turning_points.bounces = particle.orbit.turning_points.N - 1;

  % Find the indexes in the time array where these are:
  for n = 1:particle.orbit.turning_points.N
    particle.orbit.turning_points.idx(n) = max(find(particle.orbit.t <= particle.orbit.turning_points.t(n)));
  end  
  particle.orbit.turning_points.x = interp1(particle.orbit.t, particle.orbit.position.x, particle.orbit.turning_points.t);
  particle.orbit.turning_points.y = interp1(particle.orbit.t, particle.orbit.position.y, particle.orbit.turning_points.t);
  particle.orbit.turning_points.z = interp1(particle.orbit.t, particle.orbit.position.z, particle.orbit.turning_points.t);
  particle.orbit.turning_points.R = sqrt(particle.orbit.turning_points.x.^2 + particle.orbit.turning_points.y.^2);

  particle.orbit.turning_points.bounce_Bmax = interp1(particle.orbit.t, particle.orbit.magnetic_field.B, particle.orbit.turning_points.t);
  particle.orbit.turning_points.bounce_Bmin(1) = min(particle.orbit.magnetic_field.B(1:particle.orbit.turning_points.idx(1)));
  if (particle.orbit.turning_points.N > 1)
     for n = 2:particle.orbit.turning_points.N
         particle.orbit.turning_points.bounce_Bmin(n) = min(particle.orbit.magnetic_field.B(particle.orbit.turning_points.idx(n-1):particle.orbit.turning_points.idx(n)));
     end
  end
  particle.orbit.turning_points.mirror_ratio = particle.orbit.turning_points.bounce_Bmax./particle.orbit.turning_points.bounce_Bmin;

  % Calculate the bounce frequency
  if (particle.orbit.turning_points.bounces >= 1)

    particle.orbit.turning_points.bounce_period = diff(particle.orbit.turning_points.t);
    particle.orbit.turning_points.bounce_frequency = 1./particle.orbit.turning_points.bounce_period;
    particle.orbit.turning_points.bounce_angular_frequency = 2*pi*particle.orbit.turning_points.bounce_frequency;
  end
else
  particle.orbit.turning_points.status = 'passing';
  particle.orbit.classification.label2 = 'Not bouncing';
  particle.orbit.classification.status2 = 0;  
  
end

printf('%f s\n', toc()); fflush(stdout());

%{


function [particle] = MOrbit_turning_points(particle)
% ------------------------------------------------------------------------
% Turning points location, bounce frequency and mirror ratio
% ------------------------------------------------------------------------
counter_passing = 1;
counter_trapped = 0;

for h=0: (length(particle.guidingcenter.t)-2)
  if (particle.guidingcenter.position.phi(h+1) - particle.guidingcenter.position.phi(h+2) < 0.01)
    counter_passing = counter_passing + 1;
  elseif (particle.guidingcenter.position.phi(h+1) - particle.guidingcenter.position.phi(h+2) > 0.01)
    counter_trapped = counter_trapped + 1;
  endif
endfor
counter_passing;
counter_trapped;

if (counter_trapped >0)

  elseif (counter_trapped  == 0);

endif




%Old code




tic();
printf('%s', 'Turning points: '); fflush(stdout());
% ORBIT
% Search for where lambda = 0, that is when the parellel
% velocity is zero: lambda = cos(vpll/v)
particle.orbit.turning_points.t = zerocrossing(particle.orbit.t, particle.orbit.lambda);

% Calculate the location of the turning points
if (isempty(particle.orbit.turning_points.t) != 1)

  particle.orbit.turning_points.status = 'trapped';
  particle.orbit.classification.label2 = 'Bouncing';
  particle.orbit.classification.status2 = 1;
  particle.orbit.turning_points.N = length(particle.orbit.turning_points.t);
  particle.orbit.turning_points.bounces = particle.orbit.turning_points.N - 1;

  % Find the indexes in the time array where these are:
  for n = 1:particle.orbit.turning_points.N
    particle.orbit.turning_points.idx(n) = max(find(particle.orbit.t <= particle.orbit.turning_points.t(n)));
  end  
  particle.orbit.turning_points.x = interp1(particle.orbit.t, particle.orbit.position.x, particle.orbit.turning_points.t);
  particle.orbit.turning_points.y = interp1(particle.orbit.t, particle.orbit.position.y, particle.orbit.turning_points.t);
  particle.orbit.turning_points.z = interp1(particle.orbit.t, particle.orbit.position.z, particle.orbit.turning_points.t);
  particle.orbit.turning_points.R = sqrt(particle.orbit.turning_points.x.^2 + particle.orbit.turning_points.y.^2);

  particle.orbit.turning_points.bounce_Bmax = interp1(particle.orbit.t, particle.orbit.magnetic_field.B, particle.orbit.turning_points.t);
  particle.orbit.turning_points.bounce_Bmin(1) = min(particle.orbit.magnetic_field.B(1:particle.orbit.turning_points.idx(1)));
  if (particle.orbit.turning_points.N > 1)
     for n = 2:particle.orbit.turning_points.N
         particle.orbit.turning_points.bounce_Bmin(n) = min(particle.orbit.magnetic_field.B(particle.orbit.turning_points.idx(n-1):particle.orbit.turning_points.idx(n)));
     end
  end
  particle.orbit.turning_points.mirror_ratio = particle.orbit.turning_points.bounce_Bmax./particle.orbit.turning_points.bounce_Bmin;

  % Calculate the bounce frequency
  if (particle.orbit.turning_points.bounces >= 1)

    particle.orbit.turning_points.bounce_period = diff(particle.orbit.turning_points.t);
    particle.orbit.turning_points.bounce_frequency = 1./particle.orbit.turning_points.bounce_period;
    particle.orbit.turning_points.bounce_angular_frequency = 2*pi*particle.orbit.turning_points.bounce_frequency;
  end
else
  particle.orbit.turning_points.status = 'passing';
  particle.orbit.classification.label2 = 'Not bouncing';
  particle.orbit.classification.status2 = 0;  
  
end

printf('%f s\n', toc()); fflush(stdout());
%}