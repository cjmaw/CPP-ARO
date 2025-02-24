function [midplaneStrainsBody,midplaneStrainsPrincipal] = LaminateMidplaneStrains(eK,t,angles)
% Computes mid-thickness strains for each ply in laminate from laminate
% midplane strains and curvatures, lamina thicknesses, and transformation 
% matrices of each ply concatenated vertically.

% Author: Cade Maw
% Date: 2/23/25

fprintf('\n<strong>--- COMPUTE LAMINATE MIDPLANE STRAINS ---</strong>\n\n')

strains = eK(1:3,1);
curvatures = eK(4:6,1);
N = numel(t); % number of layers

%% Display inputs
fprintf('\n<strong>INPUTS:</strong>\n\n')

% Thicknesses
fprintf('Thicknesses (inches):\n\n')
for i = 1:N
    fprintf('{ t_%2.0f } = { %5.4f }\n',i,t(i))
end
fprintf('\n')

% Strains and Curvatures
fprintf('Strains and Curvatures:')
Pretty6x1(eK,'e','K','uin/in','urad/in')
fprintf('\n')

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

%% Compute transformation matrices for each ply
TList = zeros(3*N,3);
for i = 1:N
    [TList((3*i-2:3*i),1:3),~] = transform(angles(i));
end


%% Compute midplane strains for each ply
midplaneStrainsBody = zeros(3,N);
for i = 1:N
    midplaneStrainsBody(:,i) = strains + zbar(i).*curvatures;
end

midplaneStrainsPrincipal = zeros(3,N);
R = reuter();
for i = 1:N
    midplaneStrainsPrincipal(:,i) = R*TList((3*i-2:3*i),1:3)*inv(R)*midplaneStrainsBody(:,i);
end

%% Print midplane strains for each ply

fprintf('\n<strong>RESULTS:</strong>')
fprintf('\nMidplane Strains in Body Coordinates')
for i = 1:N
    fprintf('\n\nPly No: %2.0f',i)
    Pretty3x1(midplaneStrainsBody(:,i),'e','in/in')
end
fprintf('\n')

fprintf('\n\nMidplane Strains in Principal Coordinates')
for i = 1:N
    fprintf('\n\nPly No: %2.0f',i)
    Pretty3x1(midplaneStrainsPrincipal(:,i),'e','in/in')
end
fprintf('\n')

return