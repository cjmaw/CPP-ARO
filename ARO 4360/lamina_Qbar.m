function [Qbarmat] = lamina_Qbar(E1,E2,G12,v12,x)
% Computes reduced stiffness matrix Qbar given lamina material properties
% and angle between material principal axes and body coordinate axes.
% Angle must be in DEGREES.

% Author: Cade Maw
% Date Created: 2/4/2025
% Last Updated: 2/4/2025

% Compute T^-1, Q, R, T, R^-1
[T,Ti] = transform(x);
Q = lamina_Q(E1,E2,G12,v12);
R = reuter();
Ri = inv(R);

% Compute Qbar
Qbarmat = Ti*Q*R*T*Ri;

return
