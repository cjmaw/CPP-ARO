function CantileverShearMoment ()
% Generates Shear and Bending Moment diagrams for a cantilever beam with
% point loads. 
% IMPORTANT NOTE: DOWNWARD POINT LOAD IS POSITIVE

% Author: Cade Maw
% Date Created: 2/15/25
% Last Updated: 2/15/25

%% ---- INPUTS AND SETUP ----
% Change this stuff
    file = 'MainWing_a=0.00_v=51.22fts.xlsx';
    rangeLIFT = 'V59:V95';
    rangeSPAN = 'O59:O95';
    lift = 3.*readmatrix(file,'Range',rangeLIFT);
    span = readmatrix(file,'Range',rangeSPAN);

length = 4.411;
forces = -lift; % DOWNWARD POINT LOAD IS POSITIVE
locs = span;

% Don't change this stuff
numLoads = numel(forces);

fprintf('\n<strong>INPUTS:</strong>\n')
InputsTable = table(locs,forces,'VariableNames',{'Locations (ft)','Forces (lb)'});
disp(InputsTable)
fprintf('\nNOTE: Downward point load is positive.\n')

%% ---- SOLVE FOR REACTIONS ----
P_rxn = sum(forces);
BM_rxn = sum(forces.*locs);

fprintf('\n<strong>CANTILEVER REACTIONS:</strong>\n')
fprintf('Reaction force = %5.3f lb (+ is up)',P_rxn)
fprintf('\n')
fprintf('Reaction moment = %4.3f ft-lb (+ is clockwise)',BM_rxn)
fprintf('\n')


%% ---- COMPUTE SHEAR ----
V = zeros(1,numLoads+2);
V(1) = P_rxn;
V(2) = P_rxn - forces(1);
V(numLoads+2) = 0;
for i = 3:numLoads
    V(i) = V(i-1) - forces(i-1);
end


%% ---- COMPUTE BENDING MOMENT ----
M = zeros(1,numLoads+2);
M(1) = -BM_rxn;
M(2) = -BM_rxn + P_rxn*locs(1);
M(numLoads+1) = 0;
M(numLoads+2) = 0;

for i = 3:numLoads
    j = i-2;
    sub = zeros(1,j);
    for idx = 1:j
        sub(idx) = forces(idx)*(locs(i-1)-locs(idx));
    end
    totalSubtract = sum(sub);
    M(i) = -BM_rxn + P_rxn*(locs(i-1)) - totalSubtract;
end


%% ---- PLOTS ----

xplot = [locs;length];

% Plot shear distribution
figure(6); hold on
line([0 xplot(1)],[V(1) V(1)],'Color','k')
    fill([0 xplot(1) xplot(1) 0],[0 0 V(1) V(1)],'b','FaceAlpha',0.25,'EdgeColor','none')
line([xplot(1) xplot(1)],[V(1) V(2)],'Color','k')
    for i = 2:numLoads
        line([xplot(i-1) xplot(i)],[V(i) V(i)],'Color','k')
            fill([xplot(i-1) xplot(i) xplot(i) xplot(i-1)],[0 0 V(i) V(i)],'b','FaceAlpha',0.25,'EdgeColor','none')
        line([xplot(i) xplot(i)],[V(i+1),V(i)],'Color','k')
    end
line([xplot(numLoads) xplot(numLoads+1)],[V(end) V(end)],'Color','k')
title('Shear Force Diagram')
xlabel('Location (ft)')
ylabel('Shear (lb)')
axis padded
hold off
grid on

xplot = [0;xplot];

% Plot bending moment distribution
figure(7); hold on
plot(xplot,M,'Color','k')
area(xplot,M,'EdgeColor','k','FaceColor','b','FaceAlpha',0.25)
title('Bending Moment Diagram')
xlabel('Location (ft)')
ylabel('Bending Moment (ft-lb)')
hold off
grid on

return