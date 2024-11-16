clear; clc

mu = 398600.435; % km^3/s^2 for Earth

radius_LEO = 6928; % km
radius_GEO = 42378; % km
a_transfer = (radius_LEO+radius_GEO)/2;

equation_1 = (2/radius_LEO) - (1/a_transfer);
velocity_transfer = sqrt(mu*equation_1);

TXfer = 2*pi*sqrt(a_transfer^3/mu); % s
dt = round(TXfer/100);

q0 = [radius_LEO;0;0;0;velocity_transfer;0];
q(:,1)=q0;

for i = 1:62 % Can pick different values here to get it all the way around the orbit
    f(i) = 1 - ( mu*dt^2 / (2*norm(q(1:3,i))^3) ) + (mu*dt^3 / (2*norm(q(1:3,i))^5 )) * dot(q(1:3,i),q(4:6,i));
    g(i) = dt - (mu*dt^3 / (6*norm(q(1:3,i))^3) );
    f_dot(i) = -mu*dt/(norm(q(1:3,i))^3) + (3*mu*dt^2/(2*norm(q(1:3,i))^5)) * dot(q(1:3,i),q(4:6,i));
    g_dot(i) = 1 - mu*dt^2/(2*norm(q(1:3,i))^3);

    q(1:3,i+1) = f(i).*q(1:3,i) + g(i).*q(4:6,i);
    q(4:6,i+1) = f_dot(i).*q(1:3,i) + g_dot(i).*q(4:6,i);
end

lagrangeCoeff_x = q(1,:);
lagrangeCoeff_y = q(2,:);

figure
plot(lagrangeCoeff_x,lagrangeCoeff_y)