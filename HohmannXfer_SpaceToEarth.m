clear; clc

mu = 398600.435; % for Earth

%   Coming from deep space at A
V_A = 10; % km/s
R_A = 5000+6378; % km

%   Arriving at 500 km altitude at B
R_B = 500+6378; % km
V_B = sqrt(mu/R_B); % km/s

%   Transfer orbit
a_Xfer = (R_A+R_B)/2;
V_Xfer_A = sqrt(mu*( (2/R_A)-(1/a_Xfer) ));
V_Xfer_B = sqrt(mu*( (2/R_B)-(1/a_Xfer) ));

deltaV_A = abs(V_Xfer_A-V_A);
deltaV_B = abs(V_Xfer_B-V_B);

TotalDeltaV = deltaV_B+deltaV_A;
fprintf('Total DV is %f km/s.',TotalDeltaV)