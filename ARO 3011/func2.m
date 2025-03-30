function f2 = func2(t,u1,u2)
    gamma = 1.4;
    Vr = u1;
    Vt = u2;

    val1 = Vt^2*Vr - (((gamma-1)/2)*(1-Vr^2-Vt^2)*(2*Vr+Vt*cot(t)));
    val2 = ((gamma-1)/2)*(1-Vr^2-Vt^2) - Vt^2;
    
    f2 = val1/val2;
return