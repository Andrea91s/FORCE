function [gradB_x, gradB_y, gradB_z, gradBmod] = gradBnew(x,y,z)
% ------------------------------------------------------------------------
% gradB calculation
% ------------------------------------------------------------------------

  % Calculate  the gradient box
  ds = 0.01;	% grid step in X, Y and Z in m
  NB = 2;
  s = linspace(-NB*ds, NB*ds, 2*NB+1);
  N = length(x);
  M = length(s);
  
  % Convert the array s into a matrix M x N
  sm = repmat(s', 1, N);
  
  % Convert the array x,y,z into a matrix M x N
  xm = repmat(x, M, 1);
  ym = repmat(y, M, 1);
  zm = repmat(z, M, 1);
  
  % Generate a matrix in which the rows are the individual x(i) to which
  % the s have been added:
  SX = sm + xm;
  SY = sm + ym;
  SZ = sm + zm;
  
  % Now evaluate the field along the columns
  for m = 1:M
     [BX, BY, BZ, BBX(m,:)] = BfieldFast(SX(m,:), y, z);
     [BX, BY, BZ, BBY(m,:)] = BfieldFast(x, SY(m,:), z);
     [BX, BY, BZ, BBZ(m,:)] = BfieldFast(x, y, SZ(m,:));
  end
  
  % And now evaluate the gradient along s
  for n = 1:N
    GB_x = gradient(BBX(:,n), ds); 
    GB_y = gradient(BBY(:,n), ds); 
    GB_z = gradient(BBZ(:,n), ds); 
    gradB_x(n) = GB_x(NB+1);
    gradB_y(n) = GB_y(NB+1);
    gradB_z(n) = GB_z(NB+1);
  end
  gradBmod = sqrt(gradB_x.^2 + gradB_y.^2 + gradB_z.^2);   
endfunction
