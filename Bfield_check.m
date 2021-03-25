% Script to test the toroidal symmetry of the B field

R = linspace(0.1, 2, 100);
Z = 0.5;
PHI1 = 0; X1 = R*cos(PHI1); Y1 = R*sin(PHI1); Z1 = Z;
PHI2 = pi/2; X2 = R*cos(PHI2); Y2 = R*sin(PHI2); Z2 = Z;
PHI3 = pi; X3 = R*cos(PHI3); Y3 = R*sin(PHI3); Z3 = Z;
PHI4 = 3*pi/2; X4 = R*cos(PHI4); Y4 = R*sin(PHI4); Z4 = Z;
R1 = sqrt(X1.^2 + Y1.^2);
R2 = sqrt(X2.^2 + Y2.^2);
R3 = sqrt(X3.^2 + Y3.^2);
R4 = sqrt(X4.^2 + Y4.^2);

[BX1, BY1, BZ1, BB1] = BfieldFast(X1, Y1, Z1);
[BX2, BY2, BZ2, BB2] = BfieldFast(X2, Y2, Z2);
[BX3, BY3, BZ3, BB3] = BfieldFast(X3, Y3, Z3);
[BX4, BY4, BZ4, BB4] = BfieldFast(X4, Y4, Z4);

figure(1, 'position', [100 100 1000 460]); 
subplot(1,4,1)
plot(R1,BX1,R2,BX2,R3,BX3,R4,BX4)
xlabel('major radius R (m)', 'fontsize', 12);
ylabel('BX (Tesla)')
legend('0', '\pi/2', '\pi', '3 \pi/2')
title(['Z = ' num2str(Z) '(m)'])

subplot(1,4,2)
plot(R1,BY1,R2,BY2,R3,BY3,R4,BY4)
xlabel('major radius R (m)', 'fontsize', 12);
ylabel('BY (Tesla)')
legend('0', '\pi/2', '\pi', '3 \pi/2')
title(['Z = ' num2str(Z) '(m)'])

subplot(1,4,3)
plot(R1,BZ1,R2,BZ2,R3,BZ3,R4,BZ4)
xlabel('major radius R (m)', 'fontsize', 12);
ylabel('BZ (Tesla)')
legend('0', '\pi/2', '\pi', '3 \pi/2')
title(['Z = ' num2str(Z) '(m)'])

subplot(1,4,4)
plot(R1,BB1,R2,BB2,R3,BB3,R4,BB4)
xlabel('major radius R (m)', 'fontsize', 12);
ylabel('B TOTAL (Tesla)')
legend('0', '\pi/2', '\pi', '3 \pi/2')
title(['Z = ' num2str(Z) '(m)'])

