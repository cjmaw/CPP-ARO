function [QbarAssembly] = AssembleQbars(E1,E2,G12,v12,angles)
% This program computes Qbar matrices for all plies and concatenates them
% into a 3*N x 3 matrix for use in ABBD.m program.
% All inputs should be column vectors, with the top ply at index 1.

% Author: Cade Maw
% Date: 2/24/25

N = numel(E1);

            if numel(E2) ~= N
                fprintf('Error in number of elements in E2')
            end
            if numel(G12) ~= N
                fprintf('Error in number of elements in G12')
            end
            if numel(v12) ~= N
                fprintf('Error in number of elements in v12')
            end
            if numel(angles) ~= N
                fprintf('Error in number of elements in angles')
            end

QbarAssembly = zeros(3*N,3);
for i = 1:N
    QbarAssembly(3*i-2:3*i,1:3) = Qbar(E1(i),E2(i),G12(i),v12(i),angles(i));
end

return