% Homework 3
clear; clc

%% Inputs
M = 0.7;
xbarcg = 0.3;
alphaZL = -3.*pi/180; % rad
Cmac = -0.01;
etah = 0.9;
etav = 

% Control Surface Effectiveness Factors
taue = 0.44; % Elevator
taur = ; % Rudder
taua = ; % Aileron

SF = 196/2.4; % Scale factor from paper measurements

% Wing and Tail Distances
xwh = SF*1.6; % Longitudinal distance from wing LE to HT LE
ywh = SF*0.15; % Vertical distance from top of wing to top of HT

% Wing
bw = 196;
crw = SF*0.64;
ctw = SF*0.18;
xw = SF*1.05;
taperw = ctw/crw;
sw = (bw/2)*crw*(1+taperw);
ARw = bw^2 / sw;
cbarw = (2/3)*crw*(1+taperw+taperw^2)/(1+taperw);
LEanglew = atan(xw/(bw/2)); % rad
xbarac_w = 0.25;
dihedralw = % radians
ymgcw = 

% Horizontal Tail
bh = SF*0.89;
crh = SF*0.4;
cth = SF*0.1;
xh = SF*0.4;
taperh = cth/crh;
sh = (bh/2)*crh*(1+taperh);
ARh = bh^2 / sh;
cbarh = (2/3)*crh*(1+taperh+taperh^2)/(1+taperh);
LEangleh = atan(xh/(bh/2)); % rad
xmgc_h = bh*(1+2*taperh)*tan(LEangleh) / (6*(1+taperh));
xmgc_w = bw*(1+2*taperw)*tan(LEanglew) / (6*(1+taperw));
xbarac_h = (1/cbarw)*(xwh - xmgc_w + xmgc_h + 0.25*cbarh);

% Vertical Tail
bv = 
crv = 
ctv = 
xv = 
taperv = 
sv = (bv/2)*crv*(1+taperv);
ARv = bv^2 / sv;
ARveff = 2*ARv;
cbarv = (2/3)*crv*(1+taperv+taperv^2)/(1+taperv);
LEanglev = atan(xh/(bv/2)); % rad
ymcgv = 

% Other 
r = xwh/(bw/2);
m = ywh/(bw/2);

zw = % IDK WHAT THIS IS

%% LONGITUDINAL DIRECTION
% Lift Coeff
CLaw = polhamus(ARw,taperw,LEanglew,M);
CLo = CLaw*abs(alphaZL);
deda = downwash(ARw,taperw,LEanglew,M,r,m);
CLah = polhamus(ARh,taperh,LEangleh,M);
CLa_AC = CLaw + etah*(sh/sw)*CLah*(1-deda);
CLih = etah*(sh/sw)*CLah;
CLde = etah*(sh/sw)*CLah*taue;

% Moment Coeff
Cmo = Cmac + CLo*(xbarcg - xbarac_w);
Cma = CLaw*(xbarcg - xbarac_w) - CLah*etah*(sh/sw)*(1-deda)*(xbarac_h-xbarcg);
Cmde = -CLah*etah*(sh/sw)*taue*(xbarac_h-xbarcg);
Cmih = -CLah*etah*(sh/sw)*(xbarac_h-xbarcg);

% Neutral point
num = xbarac_w + (CLah/CLaw)*etah*(sh/sw)*(1-deda)*xbarac_h;
den = 1 + (CLah/CLaw)*etah*(sh/sw)*(1-deda);
xbarnp = num/den;

% Static Margin
SM = xbarnp - xbarcg;


% Print Results
fprintf('+--------------------------------+\n')
fprintf('|     LONGITUDINAL DIRECTION     |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('| Parameter |    Value   | Units |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('|    CLo    | %10.4f |       |\n',CLo)
fprintf('+-----------+------------+-------+\n')
fprintf('|    CL \x3B1   | %10.4f | 1/rad |\n',CLa_AC)
fprintf('+-----------+------------+-------+\n')
fprintf('|    CLih   | %10.4f | 1/rad |\n',CLih)
fprintf('+-----------+------------+-------+\n')
fprintf('|    CLde   | %10.4f | 1/rad |\n',CLde)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cmo    | %10.4f |       |\n',Cmo)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cm \x3B1   | %10.4f | 1/rad |\n',Cma)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cmde   | %10.4f | 1/rad |\n',Cmde)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cmih   | %10.4f | 1/rad |\n',Cmih)
fprintf('+-----------+------------+-------+\n\n')


fprintf('<strong>Nondimensional Neutral Point Location:</strong> %5.3f \n',xbarnp)
fprintf('\n<strong>Longitudinal Static Margin:</strong> %5.3f \n',SM)

%% LATERAL DIRECTION

lambda025 = sweepx(LEanglew,taperw,0.25);
nv1dadb = 0.724 + 3.06*(sv/sw)/(1+cos(lambda025)) + 0.4*zw/d + 0.009*ARw;
CLav = polhamus(ARveff,taperv,LEanglev,M);
Cybv = nv1dadb*(sv/sw)*CLav;
Cyb = 1.3*Cybv;

Cydr = etav*(sv/sw)*CLav*taur;

Cybw = -0.0001*dihedralw*57.3;

vv = (sv/sw)*(ymgcv/bw);
Clbv = -vv*etav*nv1dadb;
Clbw = -CLaw*dihedralw*(ymgcw/bw);
Clb = Clbv + Clbw;
Clda ...

Cldr = CLav*etav*(sv/sw)*taur*ymgcv;

Cnb = etav*(sv/sw)*nv1dadb*CLav*(xacv-xcg)/bw;
Cndr = -etav*(sv/sw)*CLav*taur*(xacv-xcg)/bw;
