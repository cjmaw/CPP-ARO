function [Rmat] = reuter()
% Writes Reuter's Matrix.
% Used to convert between true strain and engineering strain.

% Author: Cade Maw
% Date Created: 2/4/2025
% Last Updated: 2/4/2025

Rmat = [1,0,0;...
        0,1,0;...
        0,0,2];

return