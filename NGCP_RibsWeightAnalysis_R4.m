clear; clc

% ------- ARCH 1 -------
totalSpan = 14; % ft
untaperedSpan = 10; % ft
cRoot = 2.1; % ft 
cTip = 0.656; % ft
tw = 0.165; % ft ---- wing thickness
tr = 0.25/12; % ft ---- rib thickness

rpf = 2; % ribs/ft

%       Calculate Untapered Ribs Volume
numRibs_untapered = round(rpf*untaperedSpan); % round to a whole number
volume_untapered = cRoot*tw*tr*numRibs_untapered;

%       Calculate Tapered Ribs Volume
numRibs_tapered = round(totalSpan*rpf) - numRibs_untapered;

for i = 1:numRibs_tapered/2+1
    
    x(i) = (i-1)/rpf;
    c(i) = cRoot - 2*x(i)*(cRoot-cTip) / ((totalSpan-untaperedSpan)/2);
    volume(i) = c(i)*tr*tw;

end

format longg
totalVolumeArch1 = 2*(sum(volume)+volume_untapered)


% ===============================================


clear

% ------- ARCH 3 ------- 
span = 13.5; % ft
cRoot = 1.748; % ft 
cTip = 0.874; % ft
tw = 0.12*(cRoot+cTip)/2; % ft ---- wing thickness
    % ^^ Average of root and tip thicknesses for 12% thickness ratio
tr = (1/8)/12; % ft ---- rib thickness
sweepLE = 0.095*pi/180; % RADIANS

h = cRoot - cTip - span*tan(sweepLE)/2;
sweepTE = atan(2*h/span); % RADIANS

rpf = 2; % ribs/ft
numRibs = round(rpf*span)-1; % round to even #

for i = 1:numRibs/2+1
    
    x(i) = (i-1)/rpf;
    c(i) = cRoot - x(i)*tan(sweepLE) - x(i)*tan(sweepTE);
    volume(i) = c(i)*tr*tw;

end

format longg
totalVolumeArch3 = 2*sum(volume) % x2 is for both halves of wing



% ===============================================


clear

% ------- ARCH 1 with tapered trapezoidal wing ------- 
span = 14.4; % ft
cRoot = 2.176; % ft 
cTip = 0.653; % ft
tw = 0.12*(cRoot+cTip)/2; % ft ---- wing thickness
    % ^^ Average of root and tip thicknesses for 12% thickness ratio
tr = (1/8)/12; % ft ---- rib thickness
sweepLE = 3*pi/180; % RADIANS

h = cRoot - cTip - span*tan(sweepLE)/2;
sweepTE = atan(2*h/span); % RADIANS

rpf = 2; % ribs/ft
numRibs = round(rpf*span)-1; % round to even #

for i = 1:numRibs/2+1

    x(i) = (i-1)/rpf;
    c(i) = cRoot - x(i)*tan(sweepLE) - x(i)*tan(sweepTE);
    volume(i) = c(i)*tr*tw;

end

format longg
totalVolumeArch1Tapered = 2*sum(volume) % x2 is for both halves of wing