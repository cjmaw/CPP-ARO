function deda = downwash(AR,taperRatio,LEangle,M,r,m)

KA = 1/AR - 1/(1+AR^1.7);
Kl = (10-3*taperRatio)/7;
Kmr = (1-m/2)/r^0.33;

lambda25 = atan( tan(LEangle) - 4*0.25*(1-taperRatio)/(AR*(1+taperRatio)) );

val = 4.44*(KA*Kl*Kmr*sqrt(cos(lambda25)))^1.19;

deda = val/sqrt(1-M^2);
