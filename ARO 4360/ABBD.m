function ABBD (t,theta,QbarList)
% Takes column vectors for t = lamina thicknesses (inches), theta = ply
% angles (degrees), and concatenated Qbar matrices (Msi) for laminate to
% compute ABBD matrix (ksi). 
% All input vectors and matrices must be in order from top ply to bottom.

% Author: Cade Maw
% Date: 2/17/25


N = numel(t); % number of layers

% Reorganize Qbar's from a 3x(3*N) matrix to a 3x3xN matrix
QbarLaminate = zeros(3,3,N);
for k = 1:N
    QbarLaminate(1:3,1:3,k) = QbarList(k*3-2:k*3,:);
end


% Display inputs
fprintf('\n<strong>INPUTS:</strong>\n\n')

% Thicknesses
fprintf('Thicknesses (inches):\n')
for i = 1:N
    fprintf('{ t_%2.0f } = { %5.4f }\n',i,t(i))
end
fprintf('\n')

% Ply Angles
fprintf('Ply Angles (degrees):\n')
for i = 1:N
    fprintf('{ theta_%2.0f } = { %6.1f }\n',i,theta(i))
end
fprintf('\n')

% Qbar's
fprintf('Qbar Matrices for Each Ply:\n')
for i = 1:N
    symb = strcat('Qbar_',num2str(theta(i)));
    Pretty3x3(QbarList(i*3-2:i*3,:),symb,'Msi');
end


%% Compute z's and zbar's
tLam = sum(t);

z = zeros(N,1);
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
for i = 1:3
    for j = 1:3
        for k = 1:N
            tempA(i,j,k) = QbarLaminate(i,j,k)*t(k);
            tempB(i,j,k) = QbarLaminate(i,j,k)*t(k)*zbar(k);
            tempD(i,j,k) = QbarLaminate(i,j,k) * (t(k)*zbar(k)^2 + (1/12)*t(k)^3);
        end
        A(i,j) = sum(tempA(i,j,:));
        B(i,j) = sum(tempB(i,j,:));
        D(i,j) = sum(tempD(i,j,:));
    end
end


%% Combine to get ABBD
ABBDmat = [A,B;B,D].*1000; % ksi, assuming Qbars are in Msi

fprintf('\n<strong>RESULTS:</strong>\n\n')
for i = 1:N
    if i ~= floor(N/2)
        fprintf('             [ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
    elseif i == floor(N/2)
        fprintf('[ ABBD ]  =  ')
        fprintf('[ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
        fprintf(' ksi\n')
    end
end



return