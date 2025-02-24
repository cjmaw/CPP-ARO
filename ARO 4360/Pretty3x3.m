function Pretty3x3(matrix,symbol,units)
% Writes 3x3 matrix to command window in pretty way.

% Author: Cade Maw
% Date: 2/9/2025

% First row
fprintf('\n')
fprintf(1, '                  { %9.4f }',matrix(1,1))
fprintf(1, '  { %9.4f }',matrix(1,2))
fprintf(1, '  { %9.4f }',matrix(1,3))
fprintf('\n')
% Second row 
fprintf(1,['[ ',symbol,' ]'])
fprintf(1, '  = ')
fprintf(1, '   { %9.4f }',matrix(2,1))
fprintf(1, '  { %9.4f }',matrix(2,2))
fprintf(1, '  { %9.4f }',matrix(2,3))
fprintf(1,['  ',units])
fprintf('\n')
% Third row
fprintf(1, '                  { %9.4f }',matrix(3,1))
fprintf(1, '  { %9.4f }',matrix(3,2))
fprintf(1, '  { %9.4f }',matrix(3,3))
fprintf('\n')

return