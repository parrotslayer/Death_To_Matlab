% Takes in the Height of the Satellite and the Nadir Angle = FOV/2
% and returns the distance D from the satellite to the target (swath edge)
function output = Viewing_Geometry(H, lambda)
R_e = 6378137;  %meters

rho = asin(R_e/(R_e+H));
lambda0 = acos(R_e/(R_e+H));

n = atan((sin(rho)*sin(lambda))/(1-sin(rho)*cos(lambda)));

E = acos(sin(n)/sin(rho));

% a simple check to see if it works
shouldbezero = n + lambda + E - pi/2;

D = R_e*(sin(lambda)/sin(n));

output = D;
end