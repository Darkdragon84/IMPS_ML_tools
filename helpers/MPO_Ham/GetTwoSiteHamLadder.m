function [ H ] = GetTwoSiteHamLadder(d,J,J1,J2,h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[X,iY,Z] = su2gen(d);
Y = real(1i*iY);
I = eye(d);
I2 = eye(d*d);

H0 = kron(X,X) - kron(Y,Y) + kron(Z,Z);
HZ = kron(Z,I) + kron(I,Z);

H1 = kron(kron(X,I),kron(X,I)) - kron(kron(Y,I),kron(Y,I)) + kron(kron(Z,I),kron(Z,I)) + ...
     kron(kron(I,X),kron(I,X)) - kron(kron(I,Y),kron(I,Y)) + kron(kron(I,Z),kron(I,Z));
 
H2 = kron(kron(I,X),kron(X,I)) - kron(kron(I,Y),kron(Y,I)) + kron(kron(I,Z),kron(Z,I)) + ...
     kron(kron(X,I),kron(I,X)) - kron(kron(Y,I),kron(I,Y)) + kron(kron(Z,I),kron(I,Z));
 
Hon = J*(H0 - h*HZ);
H = 0.5*(kron(Hon,I2) + kron(I2,Hon)) + J1*H1 + J2*H2;
% H = kron(H0 - h*HZ,I2) + J1*H1 + J2*H2;
end

