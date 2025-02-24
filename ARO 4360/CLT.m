function [runningLoads] = CLT(ABBDmat,StrainCurve)
% Computes running loads & running moments in body coordinates from ABBD
% matrix and strains & curvatures in body coordinates.

% Author: Cade Maw
% Date: 2/18/25


% Inputs
fprintf('\n<strong>INPUTS:</strong>\n')
fprintf('ABBD w/ units: kips/in, kips, and in-kips\n')
for i = 1:6
    if i ~= 3
        fprintf('             [ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
    elseif i == 3
        fprintf('[ ABBD ]  =  ')
        fprintf('[ %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f ]\n',ABBDmat(i,1),ABBDmat(i,2),ABBDmat(i,3),ABBDmat(i,4),ABBDmat(i,5),ABBDmat(i,6))
    end
end

fprintf('\n')
fprintf('Strains and curvatures w/ units: in/in and rad/in')
fprintf(1, '\n  { e_x  }     { %8.3f }',StrainCurve(1,1))
fprintf(1, '\n  { e_y  }     { %8.3f } ',StrainCurve(2,1))
fprintf(1, '\n  { e_xy }     { %8.3f }',StrainCurve(3,1))
fprintf(1, '\n  { K_x  }     { %8.3f }',StrainCurve(4,1))
fprintf(1, '\n  { K_y  }     { %8.3f } ',StrainCurve(5,1))
fprintf(1, '\n  { K_xy }     { %8.3f }',StrainCurve(6,1))
fprintf('\n')

% Compute running loads
runningLoads = ABBDmat*StrainCurve;

% Results
fprintf('\n<strong>RESULTS:</strong>\n')
fprintf('Running loads w/ units: kip/in and in-kip/in')
fprintf(1, '\n  { N_x  }     { %8.3f }',runningLoads(1,1))
fprintf(1, '\n  { N_y  }     { %8.3f } ',runningLoads(2,1))
fprintf(1, '\n  { N_xy }     { %8.3f }',runningLoads(3,1))
fprintf(1, '\n  { M_x  }     { %8.3f }',runningLoads(4,1))
fprintf(1, '\n  { M_y  }     { %8.3f } ',runningLoads(5,1))
fprintf(1, '\n  { M_xy }     { %8.3f }',runningLoads(6,1))
fprintf('\n')

return