function Pretty6x6(matrix,symbol,units)
% Writes 6x6 matrix to command window in pretty way.

% Author: Cade Maw
% Date: 2/23/2025

fprintf('\n')
fprintf(['Units: ',units])
fprintf('\n')
for i = 1:6
    if i ~= 3
        fprintf('             [ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',matrix(i,1),matrix(i,2),matrix(i,3),matrix(i,4),matrix(i,5),matrix(i,6))
    elseif i == 3
        fprintf(['[ ',symbol,' ]  =  '])
        fprintf('[ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',matrix(i,1),matrix(i,2),matrix(i,3),matrix(i,4),matrix(i,5),matrix(i,6))
    end
end


return