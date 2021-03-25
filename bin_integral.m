
% *******************************************************************************************************************
% Function used to calculate the bin integrals
% *******************************************************************************************************************
function [BIT, BI] = bin_integral(classification_matrix, orbit_type, energy, de, lambda, dl, EN, PA, fid, check) 

% This function calculates the density of the fast ions in all the bins that have a specific orbit
% topology;
%
% INPUT
%   classification_matrix       2D matrix containing the orbit type by topology
%   orbit_type                  type of orbit topology for which the density is required (-1 to select all)
%   energy                      energy range of the centroid bin used in the calculation of the classification_matrix [eV]
%   de                          energy step [eV]
%   lambda                      lambda range of the centroid bin used in the calculation of the classification_matrix
%   dl                          lambda steo
%   EN                          energy array of the FID from TRANSP [eV]
%   PA                          pitch array of the FID from TRANSP
%   fid                         fast ion distribution from TRANSP at a given time and point in space
%   check                       keyword: 0 for skipping the check; 1 for checking every energy/lambda bin

% Search for orbits of type orbit_type in the classification_matrix
if (orbit_type > 0)
    [idx1, idx2] = find(classification_matrix == orbit_type);
else    % select all orbits (to calculate the total for example)
    [idx1, idx2] = find(classification_matrix != 0);       % select all elements
endif

if (isempty(idx1) == 1)
    BI = 0;
    BIT = 0;
    fprintf('No orbits of type %d.\n', orbit_type);
    return
endif

N = numel(idx1);

if (check == 1)
    W = fid;        % copy the FID in a temp variable used for checking
endif
    
% Iterate of i and j to calculate all the bin densities
for n = 1:N
    
    % Find the indexes in the FID EN and PA within the bin
    EN1(n) = min(find(EN >= energy(idx1(n)) - de/2));      % index of min value of EN in the bin
    EN2(n) = max(find(EN <= energy(idx1(n)) + de/2));      % index of max value of EN in the bin
    ei = EN1(n):EN2(n);                                 % indexes of EN within the bin
    
    PA1(n) = min(find(PA >= lambda(idx2(n)) - dl/2));      % index of min value of PA in the bin
    PA2(n) = max(find(PA <= lambda(idx2(n)) + dl/2));      % index of max value of PA in the bin
    pi = PA1(n):PA2(n);                                 % indexes of PA within the bin
    
    
    % Pick the FID within these indexes
    fid_ini_bin = fid(ei,pi);
    BI(n) = sum(fid_ini_bin(:))*de*dl;

    % Checking that I have selected the right region
    if (check == 1)
    
        W(ei,pi) = NaN;
        figure(4); 
        pcolor(PA, EN, W); shading flat; colormap jet; caxis([0 max(fid(:))]);
        title(sprintf('i = %d, j = %d, E = %f, lambda = %f', idx1(n), idx2(n), energy(idx1(n)), lambda(idx2(n))))
        
        fprintf('Given lambda = %.2f and E = %.2f keV:\n', lambda(idx2(n)), energy(idx1(n))/1000);
        fprintf('PA min = %.2f and PA max = %.2f:\n', PA(PA1(n)), PA(PA2(n)));
        fprintf('E min  = %.2f and E max  = %.2f keV:\n', EN(EN1(n))/1000, EN(EN2(n))/1000);
        fprintf('FID integral = %1.3g\n', BI(n))
        input ('ok?')
   
    endif
   
    clear ei pi fid_ini_bin;
   
end

% Calculate the total density of fast ions that have ab orbit topology of the chosen type
BIT = sum(BI);

endfunction