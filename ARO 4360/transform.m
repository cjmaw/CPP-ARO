function [Tmat,Timat] = transform(x)
% Writes transformation and inverse transformation matrices for given 
% angle x (theta in Coburn CSH). Angle must be in DEGREES

% Author: Cade Maw
% Date Created: 2/4/2025
% Last Updated: 2/4/2025

% Program three rows of T matrix separately
firstRow = [ (cosd(x))^2 , (sind(x))^2 , 2*sind(x)*cosd(x) ];
secondRow = [ (sind(x))^2 , (cosd(x))^2 , -2*sind(x)*cosd(x) ];
thirdRow = [ -sind(x)*cosd(x) , sind(x)*cosd(x) , (cosd(x))^2 - (sind(x))^2 ];

% Write Tmat and Ti mat
Tmat = [firstRow ; secondRow ; thirdRow];
Timat = inv(Tmat);

return
