function bias = computeJBias(m1,m2)

var1 = psi(1,m1)';
var11 = psi(1,m1-1)';
var2 = psi(1,m2)';
var21 = psi(1,m2-1)';

bias = (1./(m1+m2)).*(m1.^2+m2.^2 - m1.*(m1-1).*sqrt((var1+var2)./(var11 + var2)) - m2.*(m2-1).*sqrt((var1 + var2)./(var1 + var21)));