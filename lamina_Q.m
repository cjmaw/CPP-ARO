function [Q] = lamina_Q(E1,E2,G12,v12)
% Computes stiffness matrix from lamina material properties.
% Units of E1,E2,G12 must be the same.

% Author: Cade Maw
% Date Created: 2/1/25
% Last Updated: 2/1/25

% Compute secondary in-plane Poisson's ratio
v21 = v12*(E2/E1);              

% Compute stiffness matrix elements
Q11 = E1/(1-v12*v21);
Q12 = v12*E2/(1-v12*v21);
Q22 = E2/(1-v12*v21);
Q66 = G12;

% Compute stiffness matrix
Q = [Q11,Q12,0;...
     Q12,Q22,0;...
     0,  0,  Q66];

return