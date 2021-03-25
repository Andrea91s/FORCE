function [gradxF, gradyF, gradzF, gradmodF] = grad(x,y,z)
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
  
  f = @(x,y,z) x.^2 + y.^2 + z.^2;
  
  % Now evaluate the function along the columns
  for m = 1:M
     fx(m,:) = f(SX(m,:), y, z);
     fy(m,:) = f(x, SY(m,:), z);
     fz(m,:) = f(x, y, SZ(m,:));
  end
  
  % And now evaluate the gradient along s
  for n = 1:N
    dfx = gradient(fx(:,n), ds); 
    dfy = gradient(fy(:,n), ds); 
    dfz = gradient(fz(:,n), ds); 
    gradxF(n) = dfx(NB+1);
    gradyF(n) = dfy(NB+1);
    gradzF(n) = dfz(NB+1);
  end
  gradmodF = sqrt(gradxF.^2 + gradyF.^2 + gradzF.^2);   
endfunction
