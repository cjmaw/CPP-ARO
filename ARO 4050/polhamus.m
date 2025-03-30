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
