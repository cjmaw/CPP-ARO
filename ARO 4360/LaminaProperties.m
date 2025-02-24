function [E1,E2,G12,v12,v21,F1t,F2t,F1cu,F1cE,F1cS,F2c,F12] = LaminaProperties...
                (Vf,Ef,Gf,vf,Fft,Ffc,Ffs,Vm,Em,Gm,vm,Fmt,Fmc,Fms)
% Predicts lamina micromechanical stiffness and strength properties.
% Based on the Rule of Mixtures (NOT Ekvall's elastic moduli).

% Author: Cade Maw
% Date: 2/15/25


% Compute elastic moduli E1 and E2
E1 = Ef*Vf + Em*Vm;
E2_denom = Vm*Ef + Vf*Em;
E2 = Ef*Em/E2_denom;

% Compute shear moduli G12
G12_denom = Vm*Gf + Vf*Gm;
G12 = Gm*Gf/G12_denom;

% Compute Poisson's ratios
v12 = Vm*vm + Vf*vf;
v21 = v12*(E2/E1);

% Compute tensile strengths
F1t = Fft*Vf + (Em/Ef)*Fft*Vm;
F2t = Fmt;

% Compute compressive strengths
F1cu = Ffc*Vf + (Em/Ef)*Ffc*Vm;
F1cE = 2*(Vf + Vm*(Em/Ef))*sqrt(Vf*Em*Ef/(3*Vm));
F1cS = Gm/Vm;
F2c = Fmc;

% Compute shear strength
F12 = Ffs*Vf + Fms*Vm;


% Print inputs
fprintf('\n<strong>INPUTS</strong>:\n')
LaminaProperties = {'Volume Fraction';'Longitudinal Modulus (Msi)';...
                       'Shear Modulus (Msi)';'Poissons Ratio';...
                       'Tensile Strength (ksi)';'Compressive Strength (ksi)';...
                       'Shear Strength (ksi)'};
Fiber = [Vf;Ef;Gf;vf;Fft;Ffc;Ffs];
Matrix = [Vm;Em;Gm;vm;Fmt;Fmc;Fms];
InputsTable = table(LaminaProperties,Fiber,Matrix);
disp(InputsTable)

% Print outputs
fprintf('\n<strong>RESULTS</strong>:\n')
fprintf('E1 = %8.5f Msi\nE2 = %8.5f Msi\nG12 = %8.5f Msi\n',E1,E2,G12)
fprintf('v12 = %5.4f \nv21 = %8.7f \n',v12,v21)
fprintf('F1t = %6.3f ksi\nF2t = %6.3f ksi\n',F1t,F2t)
fprintf('F1cu = %6.3f ksi\nF1cE = %7.4f Msi\n',F1cu,F1cE)
fprintf('F1cS = %6.3f Msi\nF2c = %7.4f ksi\n',F1cS,F2c)
fprintf('F12 = %7.4f ksi\n',F12)