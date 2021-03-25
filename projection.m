function [ParComp, PerComp, Angle] = projection(V1, V2);
% Function that calculates the parallel and 
% perpendicular components of V1 along V2

% size of vector V1
[M, N] = size(V1);

for m = 1:M
  ParComp(m,:) = (dot(V1(m,:), V2(m,:))/norm(V2(m,:))^2)*V2(m,:);
  PerComp(m,:) = V1(m,:) - ParComp(m,:);
  Angle(m,:) = acos(dot(V1(m,:), V2(m,:))/(norm(V1(m,:))*norm(V2(m,:))));
end