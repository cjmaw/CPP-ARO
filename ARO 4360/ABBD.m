function [ABBDmat] = ABBD2 (t,theta,QbarList)
% Takes column vectors for t = lamina thicknesses (inches), theta = ply
% angles (degrees), and concatenated Qbar matrices (Msi) for laminate to
% compute ABBD matrix (ksi). 
% All input vectors and matrices must be in order from top ply to bottom.

% Author: Cade Maw
% Date: 2/17/25

fprintf('\n<strong>--- COMPUTE ABBD MATRIX ---</strong>\n\n')

N = numel(t); % number of layers


% % Display inputs
% fprintf('\n<strong>INPUTS:</strong>\n\n')
% 
% % Thicknesses
% fprintf('Thicknesses (inches):\n')
% for i = 1:N
%     fprintf('{ t_%2.0f } = { %5.4f }\n',i,t(i))
% end
% fprintf('\n')
% 
% % Ply Angles
% fprintf('Ply Angles (degrees):\n')
% for i = 1:N
%     fprintf('{ theta_%2.0f } = { %6.1f }\n',i,theta(i))
% end
% fprintf('\n')
% 
% % Qbar's
% fprintf('Qbar Matrices for Each Ply:\n')
% for i = 1:N
%     symb = strcat('Qbar_',num2str(theta(i)));
%     Pretty3x3(QbarList(i*3-2:i*3,:),symb,'Msi');
% end


%% Compute z's and zbar's
tLam = sum(t);

z = zeros(N+1,1);
z(1) = -tLam/2; % z0 in CSHv1
z(2) = -tLam/2 + t(1); %z1 in CSHv1, etc.
for i = 3:N+1
    z(i) = z(i-1) + t(i-1);
end


zbar = zeros(N,1);
for i = 1:N
    zbar(i) = z(i) + t(i)/2;
end


%% Compute A,B,D matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
for i = 1:N
    A = A + QbarList(3*i-2:3*i,1:3)*t(i);
    B = B + QbarList(3*i-2:3*i,1:3)*t(i)*zbar(i);
    D = D + QbarList(3*i-2:3*i,1:3)*(t(i)*zbar(i)^2 + (1/12)*t(i)^3);
end


%% Combine to get ABBD
ABBDmat = [A,B;B,D].*1000; % kips/in, kips, and in-kips, assuming Qbars are in Msi

fprintf('\n<strong>RESULTS:</strong>\n\n')
fprintf('Units: kips/in, kips, and in-kips\n')
for i = 1:6
    if i ~= 3
        fprintf('             [ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
    elseif i == 3
        fprintf('[ ABBD ]  =  ')
        fprintf('[ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
    end
end



return
