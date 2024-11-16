clear; clc
% ===============================
% ----- Givens -----
% ===============================
Re = 1E6;
nuAir = 1.46E-5;
alpha = -7*pi/180; % radians

% ===============================
% ----- Import Airfoil Data -----
% ===============================

airfoilData = open("AirfoilData.mat").dataValues;
trailingEdge = [1.02,-0.0223];
airfoilData = [airfoilData;trailingEdge];

Vo = Re*nuAir/trailingEdge(1);
%% =
% ======================================
% ----- Find Panel Boundary Points -----
% ======================================

N = 128; % Number of panels -- MUST BE AN EVEN NUMBER
TEx = trailingEdge(1); % X location of trailing edge


%   Upper Surface
    upperSurface = airfoilData(1:102,:);
    upperSurface = [upperSurface;trailingEdge];
    N_upper = N/2;

%       Assign x and y
        x = zeros(length(upperSurface)+1,1); % origin as first point
        y = zeros(length(upperSurface)+1,1);

        x(2:end,1) = upperSurface(:,1); 
        y(2:end,1) = upperSurface(:,2);
        xBoundaryUpper = linspace(0,TEx,N_upper+1);

%       Linear interpolation wrt x
        yBoundaryUpper = interp1(x,y,xBoundaryUpper);



%   Lower Surface
    lowerSurface = airfoilData(103:206,:);
    N_lower = N/2;

%       Assign x and y
        x = lowerSurface(:,1);
        y = lowerSurface(:,2);
        xBoundaryLower = linspace(0,TEx,N_lower+1);

%       Linear interpolation wrt x
        yBoundaryLower = interp1(x,y,xBoundaryLower);
%%
% ===============================
% ----- Find Control Points -----
% ===============================

%   Pre-allocate control point array sizes
    xControlUpper = zeros(1,N/2);
    yControlUpper = zeros(1,N/2);
    xControlLower = zeros(1,N/2);
    yControlLower = zeros(1,N/2);

%   Each control point is halfway between the neighboring boundary points
    for i = 1:N/2
        xControlUpper(i) = (xBoundaryUpper(i) + xBoundaryUpper(i+1))/2;
        yControlUpper(i) = (yBoundaryUpper(i) + yBoundaryUpper(i+1))/2;
    
        xControlLower(i) = (xBoundaryLower(i) + xBoundaryLower(i+1))/2;
        yControlLower(i) = (yBoundaryLower(i) + yBoundaryLower(i+1))/2;
    end
%%
% =============================================================
% ----- Re-Define Boundary and Control Points: CW from LE -----
% =============================================================

% Important to go clockwise (CW) to make each panel's normal vector 
% point out of the airfoil

%   Boundary Points
%       Upper Surface already goes CW from LE
%       Lower Surface goes CCW from LE - need to flip
        xBoundaryLower2 = fliplr(xBoundaryLower);
        yBoundaryLower2 = fliplr(yBoundaryLower);
    
%       Remove Duplicate Points at LE and TE
        xBoundaryLower2(1) = [];
        xBoundaryLower2(end) = [];
        yBoundaryLower2(1) = [];
        yBoundaryLower2(end) = [];
        
%       Concatenate w/ Upper first to go CW from LE
        xBoundary = [xBoundaryUpper';xBoundaryLower2'];
        yBoundary = [yBoundaryUpper';yBoundaryLower2'];

%   Control Points
%       Upper Control Points already go CW from LE
%       Lower Control Points go CCW from LE - need to flip
        xControlLower2 = fliplr(xControlLower);
        yControlLower2 = fliplr(yControlLower);

%       Concatenate w/ Upper first to go CW from LE
        xControl = [xControlUpper';xControlLower2'];
        yControl = [yControlUpper';yControlLower2']; 
%%
% =========================================
% ----- Find Panel Angles and Lengths -----
% =========================================

% Still going CW from LE

%   Make (x,y)|N+1 = (x,y)|1 so that point N can use it & complete loop
    xBoundary2 = xBoundary;
    xBoundary2(N+1) = xBoundary(1);
    yBoundary2 = yBoundary;
    yBoundary2(N+1) = yBoundary(1);

%   Pre-allocate length and angle array sizes
    theta = zeros(1,N);
    panelLength = zeros(1,N);

%   Loop to calculate angle of panel wrt x-axis and length of panel
    for i = 1:N
        theta(i) = atan2( (yBoundary2(i+1)-yBoundary2(i)) , (xBoundary2(i+1)-xBoundary2(i)) ); % radians
            % if theta(i) < 0     % Range of theta is 0 to 2pi
            %     theta(i) = theta(i) + 2*pi;
            % end
        panelLength(i) = sqrt( (yBoundary2(i+1)-yBoundary2(i))^2 + (xBoundary2(i+1)-xBoundary2(i))^2 ); 
    end
%%
% ================================
% ----- Find r|ij and phi|ij -----
% ================================

%   Make (x,y)|N+1 = (x,y)|1 so that point N can use it & complete loop
    xControl2 = xControl;
    xControl2(N+1) = xControl(1);
    yControl2 = yControl;
    yControl2(N+1) = yControl(1);

%   Pre-allocate matrices
    r = zeros(N,N);
    phi = zeros(N,N);

%   For loop to calculate r|ij and phi|ij
    for i = 1:N
        for j = 1:N
            r(i,j) = sqrt( (xControl2(j)-xControl2(i))^2 + (yControl2(j)-yControl2(i))^2 );
            phi(i,j) = atan2( (yControl2(j)-yControl2(i)) , (xControl2(j)-xControl2(i)) );
        end
    end
%%
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
            deltaS(j) = sqrt( (yBoundary2(j+1)-yBoundary2(j))^2 + (xBoundary2(j+1)-xBoundary2(j))^2 );
            C(i,j) = deltaS(j) * sin(theta(i)-phi(i,j)) / (2*pi*r(i,j));
            Cbar(i,j) = deltaS(j) * cos(theta(i)-phi(i,j)) / (2*pi*r(i,j));
        end
        B(i) = Vo*sin(theta(i)-alpha);
    end

q = C\B;
%%
% =================================
% ----- Compute Vt, Cp, and CL-----
% =================================

%   Pre-allocate matrix sizes
    Vt = zeros(N,1);
    Cp = zeros(N,1);

%   Equations from Handout
    for i = 1:N
        for j = 1:N
            Vt(i) = Vo*cos(theta(i)-alpha) - Cbar(i,j)*q(i);
            Cp(i) = 1 - (Vt(i)/Vo)^2;
        end
    end

%   Compute CL
    for i = 1:N/2
        CL1(i) = -Cp(i)*deltaS(i)*cos(theta(i)-alpha);
    end
    for i = N/2+1:N
        CL2(i) = Cp(i)*deltaS(i)*cos(theta(i)-alpha);
    end
    CL = sum(CL1)+sum(CL2);
%%
% figure
% hold on
% scatter(airfoilData(:,1),airfoilData(:,2),'MarkerEdgeColor','#D4D4D4')
% ylim([-.25 .25])
% xlim([-0.1 1.15])
% scatter(xBoundary,yBoundary,'MarkerEdgeColor','r')
% scatter(xControl,yControl,'MarkerEdgeColor','#77AC30')
% hold off
% legend('Airfoil','Boundary Points','Control Points')
