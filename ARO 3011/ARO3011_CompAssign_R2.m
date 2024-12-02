clear; clc

% ===============================
% ----- Givens -----
% ===============================

Re = 1E6;
nuAir = 1.46E-5;
alpha = -7*pi/180; % radians
N = 128; % Number of panels -- MUST BE AN EVEN NUMBER


% ===============================
% ----- Import Airfoil Data -----
% ===============================

airfoilData = open("AirfoilData.mat").dataValues;
leadingEdge = [0,0];
trailingEdge = [1.02,-0.0223];

Vo = Re*nuAir/trailingEdge(1);


% ======================================
% ----- Find Panel Boundary Points -----
% ======================================

TEx = trailingEdge(1); % X location of trailing edge


%   Upper Surface
    upperSurface = airfoilData(1:102,:);
    upperSurface = [leadingEdge;upperSurface;trailingEdge];
    N_upper = N/2;

%       Assign x and y
        x = zeros(length(upperSurface)+1,1); % origin as first point
        y = zeros(length(upperSurface)+1,1);

        x = upperSurface(:,1);
        y = upperSurface(:,2);
        xBoundaryUpper = linspace(0,TEx,N_upper+1);
        
%       Linear interpolation wrt x
        yBoundaryUpper = interp1(x,y,xBoundaryUpper);


%   Lower Surface
    lowerSurface = airfoilData(103:end,:);
    lowerSurface = [lowerSurface;trailingEdge];
    N_lower = N/2;

%       Assign x and y
        x = lowerSurface(:,1);
        y = lowerSurface(:,2);
        xBoundaryLower = linspace(0,TEx,N_lower+1);

%       Linear interpolation wrt x
        yBoundaryLower = interp1(x,y,xBoundaryLower);


% ==============================================================
% ----- Re-Define Boundary and Control Points: CCW from TE -----
% ==============================================================

%   Boundary Points
%       Upper Surface goes CW from LE - need to flip
        xBoundaryUpper2 = fliplr(xBoundaryUpper);
        yBoundaryUpper2 = fliplr(yBoundaryUpper);
%       Lower Surface already goes CCW from TE
    
%       Remove Duplicate Point at LE
%       Keep Duplicate Point at TE for Complete Loop
        xBoundaryUpper2(end) = [];
        yBoundaryUpper2(end) = [];
        
%       Concatenate w/ Upper first to go CCW from TE
        xBoundary = [xBoundaryUpper2';xBoundaryLower'];
        yBoundary = [yBoundaryUpper2';yBoundaryLower'];


% ===============================
% ----- Find Control Points -----
% ===============================

%   Pre-allocate control point array sizes
    xControl = zeros(N,1);
    yControl = zeros(N,1);

%   Each control point is halfway between the neighboring boundary points
    for i = 1:N
        xControl(i) = (xBoundary(i) + xBoundary(i+1))/2;
        yControl(i) = (yBoundary(i) + yBoundary(i+1))/2;
    end

% ***** FLIP
xControl = flip(xControl);
yControl = flip(yControl);
xBoundary = flip(xBoundary);
yBoundary = flip(yBoundary);

% =============================
% ----- Find Panel Angles -----
% =============================

%   Pre-allocate length and angle array sizes
    theta = zeros(1,N);
    panelLength = zeros(1,N);

%   Loop to calculate angle of panel wrt x-axis
    for i = 1:N
        theta(i) = atan2( (yBoundary(i+1)-yBoundary(i)) , (xBoundary(i+1)-xBoundary(i)) ); % radians
    end


% ==========================
% ----- Find r and phi -----
% ==========================

%   Pre-allocate matrices
    r = zeros(N,N);
    phi = zeros(N,N);

%   For loop to calculate r|ij and phi|ij
    for i = 1:N
        for j = 1:N
            r(i,j) = sqrt( (xControl(j)-xControl(i))^2 + (yControl(j)-yControl(i))^2 );
            phi(i,j) = atan2( (yControl(j)-yControl(i)) , (xControl(j)-xControl(i)) );
        end
    end


% =====================
% ----- Compute q -----
% =====================

%   Pre-allocate matrix sizes
    deltaS = zeros(N,1);
    C = ones(N,N).*0.5;
    Cbar = zeros(N,N);
    B = zeros(N,1);

%   Equations from Handout
    for i = 1:N
        for j = setdiff(1:N,i)    % j = 1:N but j =/= i
            deltaS(j) = sqrt( (yBoundary(j+1)-yBoundary(j))^2 + (xBoundary(j+1)-xBoundary(j))^2 );
            C(i,j) = deltaS(j) * sin(theta(i)-phi(i,j)) / (2*pi*r(i,j));
            Cbar(i,j) = deltaS(j) * cos(theta(i)-phi(i,j)) / (2*pi*r(i,j));
        end
        B(i) = Vo*sin(theta(i)-alpha);
    end

q = C\B;


% =================================
% ----- Compute Vt, Cp, and CL-----
% =================================

%   Pre-allocate matrix sizes
    Cp = zeros(N,1);
    stuff = zeros(N,1);
    CL1 = zeros(N,1);

%   Equations from Handout for Vt and Cp
    for i = 1:N
        stuff(i) = Vo*cos(theta(i)-alpha);
    end
    Vt = stuff - Cbar*q;

    for i = 1:N
        Cp(i) = 1 - (Vt(i)/Vo)^2;
    end

% Compute CL
Cp2 = Cp;
Cp2(1:N/2) = -Cp2(1:N/2);
    for i = 1:N
        CL1(1) = Cp2(i)*deltaS(i)*cos(theta(i)-alpha);
    end
CL = sum(CL1)



% figure
% hold on
% % plot(airfoilData(:,1),airfoilData(:,2),'MarkerEdgeColor','#D4D4D4')
% plot(upperSurface(:,1),upperSurface(:,2),'Color','k','LineWidth',1.6)
% plot(lowerSurface(:,1),lowerSurface(:,2),'Color','k','LineWidth',1.6)
% grid on
% ylim([-.3 .3])
% xlim([-0.1 1.15])
% xlabel('x/c')
% ylabel('y/c')
% title('SC (2)-0714 Airfoil with TE Extension')
% scatter(xBoundary,yBoundary,'MarkerEdgeColor','r')
% scatter(xControl,yControl,'MarkerEdgeColor','#77AC30')
% hold off
% legend('Airfoil','','Boundary Points','Control Points')

% Cp(1:N/2) = flip(Cp(1:N/2));
figure
hold on
plot(flip(xControl(1:N/2)),flip(Cp(1:N/2)),'DisplayName','Top')
plot(xControl(N/2+1:end),Cp(N/2+1:end),'DisplayName','Bottom')
legend('show')
title('Cp vs. Zeta (CCW)')
hold off
axis padded

% figure
% hold on
% plot(flip(Cp(1:N/2))) % should be flipped since we did CW from TE
% plot(Cp(N/2+1:N))

% T1 = table(xControl(N/2+1:N),Cp_top,Cp_bot);