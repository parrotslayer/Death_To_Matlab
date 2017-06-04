% Takes in the Height of the Satellite and the Nadir Angle = FOV/2
% and returns the distance D from the satellite to the target (swath edge)
function output = Viewing_Geometry(H, n)
R_e = 6378137;  %meters

rho = asin(R_e/(R_e+H));
lambda0 = acos(R_e/(R_e+H));

E = acos(sin(n)/sin(rho));

lambda = pi/2 - n - E;

D = R_e*(sin(lambda)/sin(n));

output = D;
end