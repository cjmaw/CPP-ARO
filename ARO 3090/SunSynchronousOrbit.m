clear;clc

J2 = 0.00108263; % for Earth
mu = 398600.435; % for Earth

%   Sun-synchronous orbit
r = 500+6378; %km
a = r; % for circular orbit
e = 0; % for circular orbit
nodeProcession = 1; % degree/day
Wdot = nodeProcession*(pi/180)/86400; % rad/s


eqn = -(3/2)*sqrt(mu)*J2*r^2 / ( ((1-e^2)^2)*a^(7/2) );
i = acos(Wdot/eqn)

fprintf('Orbit inclination for sun-synchronous orbit is %f degrees.', ...
    i*180/pi)
