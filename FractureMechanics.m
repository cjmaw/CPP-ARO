clear; clc
%                   Fracture Mechanics 
%       Cycles to Failure Using Numerical Calculus

b = 1.5; % in

f_max = 35; % ksi
f_min = -15; % ksi
delta_f = f_max - f_min; % ksi
K1C = 27; % ksi in^1/2
C = 2.09*10^(-8); % horrible units
m = 2.947;

N_step = 1;
i = 1; % iterations
a(i) = 0.01;
K1_current = BetaSECT(a(i),b)*f_max*sqrt(pi*a(i));

fprintf('Step       ai         Beta i      delta K1        N      delta a        K1         MS i\n')
fprintf('-----   ---------    ---------    ---------    -------   --------     --------    ------\n')

while K1_current < K1C
    
    B(i) = BetaSECT(a(i),b);
    delta_K1(i) = B(i) * delta_f * sqrt(pi*a(i));
    N(i) = N_step;
    delta_a(i) = N(i)*C*(delta_K1(i)^m);
    a(i+1) = a(i) + delta_a(i);
    K1(i+1) = B(i) * f_max * sqrt(pi*a(i+1));
    K1_current = K1(i+1);
    MS(i) = (K1C/K1(i+1))-1;
    
    fprintf('%i        %f     %f     %f     %i     %f     %f   %f \n',...
        i,a(i),B(i),delta_K1(i),N(i),delta_a(i),K1(i+1),MS(i))
    i = i+1;
    
end


fprintf ('Cycles to fast fracture = %i \n', N_step*(i-2))



function [beta] = BetaSECT (a, b)
    aob = a/b;
    
    num = 0.752 + 2.02*aob + 0.37*(1-sin(pi*aob/2))^3;
    den = cos(pi*aob/2);
    root = sqrt( 2*tan(pi*aob/2) / (pi*aob) );

    beta = num*root/den;

end



function [beta] = BetaCCT (a, b)
    aob = a/b;

    root = sqrt(sec(pi*aob));
    num1 = 1-0.1*aob^2 + 0.96*aob^4;
    beta = num1*root;

end 