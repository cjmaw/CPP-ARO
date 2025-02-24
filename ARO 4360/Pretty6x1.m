function Pretty6x1(vector,symbol1,symbol2,units1,units2)
% Writes 6x1 matrix to command window in pretty way.
% First three values print with symbol1 and units1.
% Last three values print with symbol2 and units2.

% Author: Cade Maw
% Date: 2/23/2025

fprintf('\n')
fprintf(1, '\n            { %9.3f } ',vector(1))
fprintf(1, '\n            { %9.3f } ',vector(2))
fprintf(1,['\n  [ ',symbol1,' ]  =  { %9.3f } ', units1],vector(3))

fprintf(1,['\n  [ ',symbol2,' ]  =  { %9.3f } ', units2],vector(4))
fprintf(1, '\n            { %9.3f } ',vector(5))
fprintf(1, '\n            { %9.3f } ',vector(6))
fprintf('\n')


return