function lambdax = sweepx(LEangle,taperRatio,x)

lambdax = atan( tan(LEangle) - 4*x*(1-taperRatio)/(AR*(1+taperRatio)) );

return