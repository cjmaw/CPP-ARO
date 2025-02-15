function PrettyStress(f,units)
% Writes 3x1 stress in body coordinates to command window in pretty way.

% Author: Cade Maw
% Date: 2/9/2025

fprintf(1, '\n  { f_xx }     { %8.3f }',f(1,1))
fprintf(1, ['\n  { f_yy }     { %8.3f } ', units],f(2,1))
fprintf(1, '\n  { f_xy }     { %8.3f }',f(3,1))
fprintf('\n')
