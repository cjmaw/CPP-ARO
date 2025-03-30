clear; clc

M1 = 5; % Upstream Mach number
gamma = 1.4; % Specific heat ratio of air
thetaShock = 12*pi/180; % Shock wave angle (rad)
h = -0.1*pi/180; % Shock wave angle step (rad)

%% Use theta-beta-M relation to find initial flow deflection angle
    % thetaS is shock wave angle, beta in Anderson
    % delta is flow deflection angle, theta in Anderson

    val = M1^2 * sin(thetaShock)^2 - 1;
    val2 = M1^2 * (gamma + cos(2*thetaShock)) + 2;
    
    delta = atan( 2*cot(thetaShock)*(val/val2) ); % Flow deflection angle (rad)

%% Use oblique shock relations to find M2

    Mn1 = M1*sin(thetaShock);
    
    val3 = Mn1^2 + (2/(gamma-1));
    val4 = (2*gamma/(gamma-1))*Mn1^2 - 1;
    Mn2 = sqrt( val3/val4 );

    M2 = Mn2 / sin(thetaShock-delta);

%% Compute Vr' and Vt' at the shock
    % Dropped the primes for code convenience
    % Vt is Vtheta

    V2 = ( (2/((gamma-1)*M2^2)) + 1 )^-0.5;
    Vr_shock = V2*cos(thetaShock-delta);
    Vt_shock = -V2*sin(thetaShock-delta);

%% Use 4th-order Runge-Kutta to solve Taylor-Maccoll Equation
    % Theta values to test, go to pi/180 since cone angle is unknown
    t = (thetaShock : h : 0).';

    % Pre-allocate matrix sizes
    k = zeros(4,2);
    w = zeros(2,numel(t));

    % Initial conditions at the shock: 
    % Vr_shock and Vt_shock when t = thetaShock
    w(1,1) = Vr_shock;
    w(2,1) = Vt_shock;

    % w1 = u1 = second derivative of Vr wrt theta
    % w2 = u2 = Vt

    % 4th-order Runge-Kutta integration loop
    for i = 1:numel(t)-1
        T = t(i);
        w1 = w(1,i);
        w2 = w(2,i);
    
        k11 = h*func1(T,w1,w2);
        k12 = h*func2(T,w1,w2);
    
        k21 = h*func1(T+0.5*h, w1+0.5*k11, w2+0.5*k12);
        k22 = h*func2(T+0.5*h, w1+0.5*k11, w2+0.5*k12);
    
        k31 = h*func1(T+0.5*h, w1+0.5*k21, w2+0.5*k22);
        k32 = h*func2(T+0.5*h, w1+0.5*k21, w2+0.5*k22);
    
        k41 = h*func1(T+h, w1+k31, w2+k32);
        k42 = h*func2(T+h, w1+k31, w2+k32);
    
        w(1,i+1) = w1 + (1/6)*(k11+2*k21+2*k31+k41);
        w(2,i+1) = w2 + (1/6)*(k12+2*k22+2*k32+k42);
    end
 
    % Boundary condition: 
    % Vt = 0 when t = thetaCone (on cone surface)
    
    % Find index where w2 is nearest to zero
    [~,idx] = min(abs(w(2,:)));

    % Clear data after cone is reached
    t(idx+1:end) = [];
    wnew = w;
    w = wnew(1:2,1:idx);
    
    % Compute cone semi-angle, thetaCone
    thetaCone = t(idx)*180/pi; % deg
    fprintf('thetaCone is %4.3f degrees.\n',thetaCone)
    
%% Compute Mach number at each theta value

    for i = 1:numel(t)
        V(i) = norm(w(1:2,i));
        M(i) = sqrt( (1/2)*(gamma-1)*((1/V(i)^2) - 1) );
    end
    
    MCone = M(idx);
    fprintf('Mach number at cone is %4.3f.\n',MCone)
    

%% Compute fluid property ratios at Stations 1 and 2

    % Isentropic flow relations at station 1
    To1_over_T1 = 1 + M1^2 * (gamma-1)/2; % To1 = To = const. everywhere
    po1_over_p1 = To1_over_T1 ^ (gamma/(gamma-1));
    ro1_over_r1 = To1_over_T1 ^ (1/(gamma-1)); % r = rho = density

    % Isentropic flow relations at station 2
    To2_over_T2 = 1 + M2^2 * (gamma-1)/2; % To2 = To = const. everywhere
    po2_over_p2 = To2_over_T2 ^ (gamma/(gamma-1)); % po2 = po = const. after shock
    ro2_over_r2 = To2_over_T2 ^ (1/(gamma-1));

    % Oblique shock relations between stations 1 and 2
    p2_over_p1 = 1 + (2*gamma)*(Mn1^2 - 1)/(gamma+1);
    r2_over_r1 = (gamma+1)*Mn1^2 / ((gamma-1)*(Mn1^2 + 2));
    T2_over_T1 = p2_over_p1 * r2_over_r1;

%% Compute fluid property ratios at each theta w/ isentropic eqns
    
    for i = 1:numel(t)
        To2_over_Ttheta(i) = 1 + M(i)^2 * (gamma-1)/2;
        po2_over_ptheta(i) = To2_over_Ttheta(i) ^ (gamma/(gamma-1));
        ro2_over_rtheta(i) = To2_over_Ttheta(i) ^ (1/(gamma-1));
    end

%% Compute fluid property ratios to plot
    % Pre-allocate array sizes
    p_over_p1 = ones(1,numel(t)+2);
    r_over_r1 = ones(1,numel(t)+2);
    T_over_T1 = ones(1,numel(t)+2);
    
    % At station 1
    p_over_p1(1) = 1;
    r_over_r1(1) = 1;
    T_over_T1(1) = 1; 

    % At station 2
    p_over_p1(2) = p2_over_p1;
    r_over_r1(2) = r2_over_r1;
    T_over_T1(2) = T2_over_T1; 

    % At all stations between 2 and cone
    for i = 3:numel(t)
        p_over_p1(i) = (1/po2_over_ptheta(i-2))*po2_over_p2*p2_over_p1;
        r_over_r1(i) = (1/ro2_over_rtheta(i-2))*ro2_over_r2*r2_over_r1;
        T_over_T1(i) = (1/To2_over_Ttheta(i-2))*To2_over_T2*T2_over_T1;
    end

%% Plots
t = [thetaCone;thetaCone-delta;t].*180./pi;
M = [M1,M2,M];

figure(1); hold on
plot(t,p_over_p1,'Color','b','DisplayName','p/p1')
plot(t,T_over_T1,'Color','r','DisplayName','T/T1')
plot(t,r_over_r1,'Color','g','DisplayName','r/r1')
plot(t,M,'Color','k','DisplayName','M')
hold off
xlabel ('Theta (deg)')
ylabel('Values')
legend show
