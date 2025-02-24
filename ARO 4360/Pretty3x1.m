function Pretty3x1(vector,symbol,units)
% Writes 3x1 matrix to command window in pretty way.

% Author: Cade Maw
% Date: 2/23/2025

fprintf(1, '\n            { %9.3f } ',vector(1))
fprintf(1,['\n  [ ',symbol,' ]  =  { %9.3f } ', units],vector(2))
fprintf(1, '\n            { %9.3f } ',vector(3))

return