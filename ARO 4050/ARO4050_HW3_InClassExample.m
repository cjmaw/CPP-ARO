% In Class Example
M = 0.2;
xbarcg = 0.264;
alphaZL = -2.*pi/180; % rad
Cmac = -0.05;
eta = 0.9;

taue = 0.6;

% Wing
bw = 36;
crw = 5.5;
ctw = 3.7;
taperw = ctw/crw;
sw = (bw/2)*crw*(1+taperw);
ARw = bw^2 / sw;
cbarw = (2/3)*crw*(1+taperw+taperw^2)/(1+taperw);
LEanglew = 0.044; % rad
xbaracw = 0.25;

% Tail
bh = 11.8;
crh = 5.8;
cth = 2.6;
taperh = cth/crh;
sh = (bh/2)*crh*(1+taperh);
ARh = bh^2 / sh;
cbarh = (2/3)*crh*(1+taperh+taperh^2)/(1+taperh);
LEangleh = 0.13; % rad
xbarach = 3.45;

% Other
r = 0.8;
m = 0.08;

% Lift Coeff
CLaw = polhamus(ARw,taperw,LEanglew,M);
CLo = CLaw*abs(alphaZL);
deda = downwash(ARw,taperw,LEanglew,M,r,m);
CLah = polhamus(ARh,taperh,LEangleh,M);
CLa_AC = CLaw + eta*(sh/sw)*CLah*(1-deda)
CLih = eta*(sh/sw)*CLah
CLde = eta*(sh/sw)*CLah*taue

% Moment Coeff
Cmo = Cmac + CLo*(xbarcg - xbaracw)
Cma = CLaw*(xbarcg - xbaracw) - CLah*eta*(sh/sw)*(1-deda)*(xbarach-xbarcg)
Cmde = -CLah*eta*(sh/sw)*taue*(xbarach-xbarcg)
Cmih = -CLah*eta*(sh/sw)*(xbarach-xbarcg)

% Neutral point
num = xbaracw + (CLah/CLaw)*eta*(sh/sw)*(1-deda)*xbarach;
den = 1 + (CLah/CLaw)*eta*(sh/sw)*(1-deda);
xbarnp = num/den

% Static Margin
SM = xbarnp - xbarcg
