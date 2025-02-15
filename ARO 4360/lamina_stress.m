function [stress_array] = lamina_stress(Q_matrix,strain_array)
% Computes stresses from lamina material properties and strains.
% Units of strain must be in/in.

% Converts between material principal coordinates if Q_matrix = stiffness matrix.
% Converts between body coordinates if Q_matrix = reduced stiffness matrix.

% Author: Cade Maw
% Date Created: 2/1/25
% Last Updated: 2/1/25

% Ensure that strain_array is a column vector
[~,c] = size(strain_array);
    if c > 1
        strain_array = strain_array';
    end

% Compute stresses
stress_array = Q_matrix * strain_array;

return
