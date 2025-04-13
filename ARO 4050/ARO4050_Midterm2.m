%% ARO4050 EXAM 2 - CADEN MAW
%%

clear; clc
%% INPUTS

M = 484/968;                                                                % Mach number
xbarcg = 0.3;                                                               % Nondimensional CG location
alphaZL = -3.*pi/180;                                                       % Zero lift AoA, in rad

Cmac = -0.01;                                                               % Zero lift moment about AC
CDo = 0.01;                                                                 % Zero-lift drag coefficient

W = 100000;                                                                 % Aircraft weight (lb)
hg = 35000;                                                                 % Geometric altitude (ft)
rho = 7.382*10^(-4);                                                        % Air density at altitude (slug/ft^3)
V = 484;                                                                    % Cruise velocity (ft/s)
e = 0.9;                                                                    % Oswald's Efficiency Factor

etah = 0.95;                                                                % Horizontal Tail Dynamic Pressure Ratio
etav = 0.95;                                                                % Vertical Tail Dynamic Pressure Ratio

taue = 0.5;                                                                 % Elevator Effectiveness Factor
taur = 0.5;                                                                 % Rudder Effectiveness Factor
taua = 0.35;                                                                % Aileron Effectiveness Factor

%% *MEASUREMENTS* 
% *Wing Dimensions*
bw = 113;                                                                   % Span   
crw = 23;                                                                   % Root chord
ctw = 4;                                                                    % Tip chord 
xLEw = 28;                                                                  % Longitudinal Dist: LE @ root to LE @ tip
zLEw = 9;                                                                   % Vertical Dist: LE @ root to LE @ tip
ailerony1 = 35;                                                             % Lateral Location of Inboard Side of Aileron
ailerony2 = 51;                                                             % Lateral Location of Outboard Side of Aileron

%%
% *Horizontal Tail Dimensions* 
bh = 46;                                                                    % Span
crh = 12;                                                                   % Root Chord
cth = 4;                                                                    % Tip chord
xLEh = 13;                                                                  % Longitudinal Dist: LE @ root to LE @ tip

%%
% *Vertical Tail Dimensions* 
bv =  31;                                                                   % Span
crv = 26;                                                                   % Root Chord
ctv = 6;                                                                    % Tip Chord
xLEv = 27;                                                                  % Longitudinal Dist: LE @ root to LE @ tip

%%
% *Wing and Tail Distances*
xwhr = 55;                                                                  % Longitudinal distance from wing LE to HT LE at root
xwh = xwhr - crw/4 + crh/4;                                                 % Longitudinal distance from wing 1/4 chord to HT 1/4 chord at root
xwv = 39;                                                                   % Longitudinal distance from wing LE to VT LE at root
zwh = 6;                                                                    % Vertical distance from top of wing to top of HT

%%
% *Other Dimensions*
zw = 0;                                                                     % Middle of fuselage to top of wing
r = xwh/(bw/2);
m = zwh/(bw/2);
d = 12;                                                                     % Maximum Diameter of Fuselage


%% CALCULATIONS FROM MEASUREMENTS
%%
% *Wing Geometry*
sweepw = atan(xLEw/((bw-d)/2));                                             % LE Sweep Angle, in rad
dihedralw = atan(zLEw/(bw/2));                                              % Dihedral Angle, in rad
taperw = ctw/crw;                                                           % Taper Ratio
sw = (bw/2)*crw*(1+taperw);                                                 % Surface Area
ARw = bw^2 / sw;                                                            % Aspect Ratio
cbarw = (2/3)*crw*(1+taperw+taperw^2)/(1+taperw);                           % MGC Length
xmgcw = bw*(1+2*taperw)*tan(sweepw) / (6*(1+taperw));                       % Longitudinal Location of MGC
ymgcw = bw*(1+2*taperw) / (6*(1+taperw));                                   % Lateral Location of MGC
xbarac_w = 0.25;                                                            % AC Location

%%
% *Horizontal Tail Geometry*
sweeph = atan(xLEh/(bh/2));                                                 % LE Sweep Angle, in rad
taperh = cth/crh;                                                           % Taper Ratio
sh = (bh/2)*crh*(1+taperh);                                                 % Surface  Area
ARh = bh^2 / sh;                                                            % Aspect Ratio
cbarh = (2/3)*crh*(1+taperh+taperh^2)/(1+taperh);                           % MGC Length
xmgch = bh*(1+2*taperh)*tan(sweeph) / (6*(1+taperh));                       % Longitudinal Location of MGC
ymgch = bh*(1+2*taperh) / (6*(1+taperh));                                   % Lateral Location of MGC
xbarac_h = (1/cbarw)*(xwhr - xmgcw + xmgch + 0.25*cbarh);                   % AC Location

%%
% *Vertical Tail Geometry*
sweepv = atan(xLEv/bv);                                                     % LE Sweep Angle, in rad
taperv = ctv/crv;                                                           % Taper Ratio
sv = (bv/2)*crv*(1+taperv);                                                 % Surface Area
ARv = bv^2 / sv;                                                            % Aspect Ratio
ARveff = 2*ARv;                                                             % Effective Aspect Ratio
cbarv = (2/3)*crv*(1+taperv+taperv^2)/(1+taperv);                           % MGC Length
xmgcv = 2*bv*(1+2*taperv)*tan(sweepv) /  (6*(1+taperv));                    % Longitudinal Location of MGC
ymgcv = 2*bv*(1+2*taperv) / (6*(1+taperv));                                 % Laterial Location of MGC
xbarac_v = (1/cbarw)*(xwv- xmgcw + xmgcv + 0.25*cbarv);                     % AC Location


%% LONGITUDINAL DIRECTION COEFFICIENTS
%%
% *Lift Coefficients*
    CLaw = polhamus(ARw,taperw,sweepw,M);
CLo = CLaw*abs(alphaZL);
    deda = downwash(ARw,taperw,sweepw,M,r,m);
    CLah = polhamus(ARh,taperh,sweeph,M);
CLa = CLaw + etah*(sh/sw)*CLah*(1-deda);
CLih = etah*(sh/sw)*CLah;
CLde = etah*(sh/sw)*CLah*taue;

%%
% *Pitching Moment Coefficients*
Cmo = Cmac + CLo*(xbarcg - xbarac_w);
Cma = CLaw*(xbarcg - xbarac_w) - CLah*etah*(sh/sw)*(1-deda)*(xbarac_h-xbarcg);
Cmde = -CLah*etah*(sh/sw)*taue*(xbarac_h-xbarcg);
Cmih = -CLah*etah*(sh/sw)*(xbarac_h-xbarcg);

%%
% *Neutral point*
    num = xbarac_w + (CLah/CLaw)*etah*(sh/sw)*(1-deda)*xbarac_h;
    den = 1 + (CLah/CLaw)*etah*(sh/sw)*(1-deda);
xbarnp = num/den;

%%
% *Static Margin*
SM = xbarnp - xbarcg;

%% LATERAL DIRECTION COEFFICIENTS
%%
% *Side Force Coefficients*
Cyo = 0;
    sweep25 = sweepx(sweepw,taperw,ARw,0.25);
    nv1dsdb = 0.724 + 3.06*(sv/sw)/(1+cos(sweep25)) + 0.4*zw/d + 0.009*ARw;
    CLav = polhamus(ARveff,taperv,sweepv,M);
    Cybv = -nv1dsdb*(sv/sw)*CLav;
    Cybw = -0.0001*abs(dihedralw)*57.3;
Cyb = 1.3*(Cybv+Cybw);
Cydr = etav*(sv/sw)*CLav*taur;
Cyda = 0;

%%
% *Rolling Moment Coefficients*
Clo = 0;
    Clbv = Cybv*ymgcv/bw;
    Clbw = -2*CLaw*dihedralw*(ymgcw/bw);
Clb = Clbv + Clbw;
    temp1 =  2*CLaw*taua*crw / (sw*bw);
    temp2 = (ailerony2^2 / 2) + (taperw-1)*(ailerony2^3 / 3) / (bw/2);
    temp3 = (ailerony1^2 / 2) + (taperw-1)*(ailerony2^3 / 3) / (bw/2);
Clda = temp1*(temp2-temp3);
Cldr = Cydr*ymgcv/bw;

%%
% *Yawing Moment Coefficients*
Cno = 0;
    xac_v = cbarw*xbarac_v;
    x_cg = cbarw*xbarcg;
Cnb = -Cybv*(xac_v-x_cg)/bw;
Cndr = -Cydr*(xac_v-x_cg)/bw;
Cnda = 0;

%%
% *Print Aerodynamic Coefficient Results*
fprintf('+--------------------------------+\n')
fprintf('|     LONGITUDINAL DIRECTION     |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('| Parameter |    Value   | Units |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('|    CLo    | %10.4f |       |\n',CLo)
fprintf('+-----------+------------+-------+\n')
fprintf('|    CL \x3B1   | %10.4f | 1/rad |\n',CLa)
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

fprintf('\n')
fprintf('+--------------------------------+\n')
fprintf('|        LATERAL DIRECTION       |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('| Parameter |    Value   | Units |\n')
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cyo    | %10.4f |       |\n',Cyo)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cy \x3B2   | %10.4f | 1/rad |\n',Cyb)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cy d \x3B1  | %10.4f | 1/rad |\n',Cyda)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cy d r  | %10.4f | 1/rad |\n',Cydr)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Clo    | %10.4f |       |\n',Clo)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cl \x3B2   | %10.4f | 1/rad |\n',Clb)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cl d \x3B1  | %10.4f | 1/rad |\n',Clda)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cl d r  | %10.4f | 1/rad |\n',Cldr)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cno    | %10.4f |       |\n',Cno)
fprintf('+-----------+------------+-------+\n')
fprintf('|    Cn \x3B2   | %10.4f | 1/rad |\n',Cnb)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cn d \x3B1  | %10.4f | 1/rad |\n',Cnda)
fprintf('+-----------+------------+-------+\n')
fprintf('|   Cn d r  | %10.4f | 1/rad |\n',Cndr)
fprintf('+-----------+------------+-------+\n\n')

fprintf('Nondimensional Neutral Point Location: %5.3f \n',xbarnp)
fprintf('\nLongitudinal Static Margin: %5.2f %% \n',SM*100)

%% TRIM CONDITIONS
% *Compute CLtrim assuming L = W*
q = 0.5*rho*V^2;                                                            % Dynamic Pressure (slug/ft^2)
CLtrim = W/(q*sw);                                                          % Trim lift coefficient

%%
% *Compute alpha_trim and ih angles*
A = [CLa, CLih; Cma, Cmih];
B = [CLtrim-CLo; -Cmo];
x = inv(A)*B;
alpha_trim = x(1);                                                          % Trim angle of attack
ih = x(2);                                                                  % Horizontal tail incidence angle

fprintf('\x3B1 trim = %5.3f deg \n',alpha_trim*180/pi)
fprintf ('ih = %5.2f deg \n',ih*180/pi)

%%
% *Compute trim thrust*
CD = CDo + (CLtrim^2) / (pi*ARw*e);
Drag = q*sw*CD;
Thrust_trim = Drag;

fprintf('Trim Thrust = %6.0f lb \n\n',Thrust_trim)

%% DETERMINING LONGITUDINAL STABILITY
%%
%     1. Cm beta should be < 0
%%
%     2. Cmo should be >= 0
%%
%     3. Static margin should be > 0
fprintf('         Cm \x3B1 = %6.4f 1/rad \n',Cma)
fprintf('         Cmo = %6.4f 1/rad \n',Cmo)
fprintf('         SM = %5.2f %% \n',SM*100)
%%
% *All three stability requirements are met, so the aircraft is longitudinally stable.*


%% DETERMINING LATERAL STABILITY
%%
%     1. Cn beta should be > 0
%%
%     2. Cl beta should be < 0
fprintf('         Cn \x3B2 = %6.4f 1/rad \n',Cnb)
fprintf('         Cl \x3B2 = %6.4f 1/rad \n',Clb)
%%
% *Both stability requirements are met, so the aircraft is laterally stable.*





%% FUNCTIONS

function deda = downwash(AR,taperRatio,LEangle,M,r,m)

KA = 1/AR - 1/(1+AR^1.7);
Kl = (10-3*taperRatio)/7;
Kmr = (1-m/2)/r^0.33;

lambda25 = atan( tan(LEangle) - 4*0.25*(1-taperRatio)/(AR*(1+taperRatio)) );

val = 4.44*(KA*Kl*Kmr*sqrt(cos(lambda25)))^1.19;

deda = val/sqrt(1-M^2);

end




function CLa = polhamus (AR,taperRatio,LEsweep,M)

if AR < 4
    k = 1 + AR*(1.87-0.000233*LEsweep)/100;
elseif AR >= 4
    k = 1 + ((8.2-2.3*LEsweep) - AR*(0.22-0.153*LEsweep))/100;
end

val = AR^2 * (1-M^2) / k^2;
tanLambda05 = tan(LEsweep) - (4*0.5*(1-taperRatio)/(AR*(1+taperRatio)));
val2 = 1 + (tanLambda05^2 / (1-M^2));

CLa = 2*pi*AR / (2 + sqrt(val*val2 + 4));

end




function lambdax = sweepx(LEangle,taperRatio,AR,x)

lambdax = atan( tan(LEangle) - 4*x*(1-taperRatio)/(AR*(1+taperRatio)) );

end