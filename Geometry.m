R_e = 6378137;
H = 6000;

rho = asin(R_e/(R_e+H));
lambda = acos(R_e/(R_e+H));

n = atan((sin(rho)*sin(lambda))/(1-sin(rho)*cos(lambda)));

E = acos(sin(n)/sin(rho));