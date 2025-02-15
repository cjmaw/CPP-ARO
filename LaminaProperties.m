function [E1,E2,G12,v12,v21] = LaminaProperties(Vf,Ef,Gf,vf,Vm,Em,Gm,vm)
% Computes lamina micromechanical stiffness properties.
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


% Print inputs
fprintf('\n<strong>INPUTS</strong>:\n')
StiffnessProperties = {'Volume Fraction';'Longitudinal Modulus';...
                       'Shear Modulus';'Poisson Ratio'};
Fiber = [Vf;Ef;Gf;vf];
Matrix = [Vm;Em;Gm;vm];
InputsTable = table(StiffnessProperties,Fiber,Matrix);
disp(InputsTable)

% Print outputs
fprintf('\n<strong>RESULTS</strong>:\n')
fprintf('E1 = %8.5f \nE2 = %8.5f \nG12 = %8.5f \n',E1,E2,G12)
fprintf('v12 = %5.4f \nv21 = %8.7f\n',v12,v21)