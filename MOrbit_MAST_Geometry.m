function [MG] = MOrbit_MAST_Geometry(fillyn, plotyn)
% Function that plot the MAST geometry poloidal cross section
% 
% INPUT
% fillyn 	keyword: 0 no fill, 1 filled

% MAST U Divertor
Divertor =  [0.715, 1.865; 0.715, 1.880; 0.8885, 1.880; 0.8885, 1.925; 1.5665, 1.925; 1.5665, 1.880; 1.6965, 1.880; 1.6965, 1.865; 1.6545, 1.865; 1.6545, 1.825; 0.7565, 1.825; 0.7565, 1.865; 0.715, 1.865];
Divertor_2 = [0.176, 0; 0.176, 1.707; 0.28, 1.707; 0.404, 1.612; 0.574, 1.612; 0.752, 1.728; 0.7835, 1.728; 0.7835, 1.716; 0.597, 1.559; 0.582, 1.550; 0.564, 1.547; 0.4165, 1.547; 0.29, 1.6735; 0.28, 1.6735; 0.28, 1.229; 0.196, 1.084; 0.196, 0];
Vessel_1 = [212.5, 29.9; 201.9, 29.9; 201.9, 105.5; 212.5, 105.6]/100;
Vessel_2 = [18, 0; 18, 200; 24, 220; 40, 220; 40, 230]/100;
Vessel_3 = [47.99, 230; 47.99, 220; 62.55, 220; 62.55, 230]/100;
Vessel_4 = [77.44, 230; 77.44, 220; 142.64, 220; 142.64, 230]/100;
Vessel_5 = [157.35, 230; 157.35, 220; 201.9, 220; 201.9, 145.9; 212.5, 145.9]/100;

% Geometry description: for Poloidal coil (copy from B field calculation file)
% Poloidal field coils
P1 = [0.115, -1.5; 0.115,1.5; 0.165, 1.5; 0.165, -1.5; 0.115, -1.5];
P2 = [0.435, 1.624; 0.435, 1.916; 0.565, 1.916; 0.565, 1.624; 0.435, 1.624];
P3radius = 0.063;
  P3centre(1,:) = [1.100, 1.100];
  theta = linspace(0, 2*pi, 61);
  P3x = P3radius*cos(theta) + P3centre(1,1);
  P3y = P3radius*sin(theta) + P3centre(1,2);
P4 = [1.405, 1.005; 1.405, 1.195; 1.595, 1.195; 1.595, 1.005; 1.405, 1.005];
P5 = [1.555, 0.405; 1.555, 0.595; 1.745, 0.595; 1.745, 0.405; 1.555, 0.405];
P6 = [1.413, 0.823; 1.413, 0.978; 1.471, 0.978; 1.471, 0.823; 1.413, 0.823];


% Plot
if (plotyn == 1)
hold on

	% Plot the divertor
	t1 = 1.5;
	plot(Divertor(:,1), Divertor(:,2), 'b', 'linewidth', t1);
	plot(Divertor(:,1), -Divertor(:,2), 'b', 'linewidth', t1);
	plot(Divertor_2(:,1), Divertor_2(:,2), 'b', 'linewidth', t1);
	plot(Divertor_2(:,1), -Divertor_2(:,2), 'b', 'linewidth', t1);
	
	% Plot the vessel
	t2 = 3;
	plot(Vessel_1(:,1), Vessel_1(:,2), 'k', 'linewidth', t2);
	plot(Vessel_1(:,1), -Vessel_1(:,2), 'k', 'linewidth', t2);
	plot(Vessel_2(:,1), Vessel_2(:,2), 'k', 'linewidth', t2);
	plot(Vessel_2(:,1), -Vessel_2(:,2), 'k', 'linewidth', t2);
	plot(Vessel_3(:,1), Vessel_3(:,2), 'k', 'linewidth', t2);
	plot(Vessel_3(:,1), -Vessel_3(:,2), 'k', 'linewidth', t2);
	plot(Vessel_4(:,1), Vessel_4(:,2), 'k', 'linewidth', t2);
	plot(Vessel_4(:,1), -Vessel_4(:,2), 'k', 'linewidth', t2);
	plot(Vessel_5(:,1), Vessel_5(:,2), 'k', 'linewidth', t2);
	plot(Vessel_5(:,1), -Vessel_5(:,2), 'k', 'linewidth', t2);
	
	% Plot the poloidal coils
	t3 = 1.5;
	plot(P1(:,1), P1(:,2), 'm', 'linewidth', t3);
	plot(P2(:,1), P2(:,2), 'm', 'linewidth', t3);
	plot(P2(:,1), -P2(:,2), 'm', 'linewidth', t3);
	plot(P3x, P3y, 'm', 'linewidth', t3) 
	plot(P3x, -P3y, 'm', 'linewidth', t3)
	plot(P4(:,1), P4(:,2), 'm', 'linewidth', t3);
	plot(P4(:,1), -P4(:,2), 'm', 'linewidth', t3);
	plot(P5(:,1), P5(:,2), 'm', 'linewidth', t3);
	plot(P5(:,1), -P5(:,2), 'm', 'linewidth', t3);
	plot(P6(:,1), P6(:,2), 'm', 'linewidth', t3);
	plot(P6(:,1), -P6(:,2), 'm', 'linewidth', t3);
	
	if (fillyn == 1)
	    fill(Divertor(:,1), Divertor(:,2), 'b')
	    fill(Divertor(:,1), -Divertor(:,2), 'b')
	    fill(Divertor_2(:,1), Divertor_2(:,2), 'b')
	    fill(Divertor_2(:,1), -Divertor_2(:,2), 'b')
	    fill(P1(:,1), P1(:,2), 'm')
	    fill(P2(:,1), P2(:,2), 'm')
	    fill(P2(:,1), -P2(:,2), 'm')
	    fill(P3x, P3y, 'm')
	    fill(P3x, -P3y, 'm')
	    fill(P4(:,1), P4(:,2), 'm')
	    fill(P4(:,1), -P4(:,2), 'm')
	    fill(P5(:,1), P5(:,2), 'm')
	    fill(P5(:,1), -P5(:,2), 'm')
	    fill(P6(:,1), P6(:,2), 'm')
	    fill(P6(:,1), -P6(:,2), 'm')
      endif

hold off
axis equal
%xlabel('R (cm)')
%ylabel('Z (cm)')
endif

function X = invertz(ZP)
[M,N] = size(ZP);
IM = ones(M,N);
IM(:,2) = -1;
X = ZP.*IM;
endfunction

MG.Vessel.part1 = Vessel_1;
MG.Vessel.part2 = Vessel_2;
MG.Vessel.part3 = Vessel_3;
MG.Vessel.part4 = Vessel_4;
MG.Vessel.part5 = Vessel_5;
MG.Coils.P1U = P1;
MG.Coils.P2U = P2;
MG.Coils.P4U = P4;
MG.Coils.P5U = P5;
MG.Coils.P6U = P6;
MG.Coils.P1L = invertz(P1);
MG.Coils.P2L = invertz(P2);
MG.Coils.P4L = invertz(P4);
MG.Coils.P5L = invertz(P5);
MG.Coils.P6L = invertz(P6);
MG.Divertor.part1 = Divertor;
MG.Divertor.part2 = Divertor_2;
MG.Divertor.part3 = invertz(Divertor);
MG.Divertor.part4 = invertz(Divertor_2);

endfunction
