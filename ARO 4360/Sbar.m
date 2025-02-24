function [Sbarmat] = Sbar (E1,E2,G12,v12,x)
% Computes reduced compliance matrix Sbar given lamina material properties
% and angle between material principal axes and body coordinate axes.
% Angle must be in DEGREES.

% Author: Cade Maw
% Date Created: 2/4/2025
% Last Updated: 2/4/2025

% Compute compliance matrix
Q = Q(E1,E2,G12,v12);
S = inv(Q);

% Compute transformation matrix for angle x
[T,~] = transform(x);

% Compute Sbar matrix
Sbarmat = (T')*S*T;

return
