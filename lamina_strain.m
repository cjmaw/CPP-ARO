function [strain_array] = lamina_strain(Q_matrix,stress_array)
% Computes strains from lamina material properties and stresses.
% Units of stress must be the same as units of E1,E2,G12.

% Converts between material principal coordinates if Q_matrix = stiffness matrix.
% Converts between body coordinates if Q_matrix = reduced stiffness matrix.

% Author: Cade Maw
% Date Created: 2/1/25
% Last Updated: 2/1/25

% Ensure that stress_array is a column vector
[~,c] = size(stress_array);
    if c > 1
        stress_array = stress_array';
    end

% Compute strains
strain_array = Q_matrix \ stress_array;

return