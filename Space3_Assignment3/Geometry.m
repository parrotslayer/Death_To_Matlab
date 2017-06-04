function output = Viewing_Geometry(H)
R_e = 6378137;  %meters
H = 6000;   %altitude of the satellite meters

rho = asin(R_e/(R_e+H));
lambda = acos(R_e/(R_e+H));

n = atan((sin(rho)*sin(lambda))/(1-sin(rho)*cos(lambda)));

E = acos(sin(n)/sin(rho));

% a simple check to see if it works
shouldbezero = n + lambda + E - pi/2;

D = R_e*(sin(lambda)/sin(n));

output = D;
end