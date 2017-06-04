% Calculates required FOV for a given Swath width and Height
% Constants
R_e = 6378137;  %m
H = 7.8335e5;   %m
Swath_Width = 200e3;  %m
deg2rad = pi/180;

% Calc variables
n = atan((Swath_Width/2)/H);

rho = asin(R_e/(R_e+H));

lambda0 = acos(R_e/(R_e+H));

E = acos(sin(n)/sin(rho));

lambda = pi/2 - n - E;

D = R_e*(sin(lambda)/sin(n));


