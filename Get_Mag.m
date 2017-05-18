% Get_Mag
% Takes in the time since epoch and the current satellite position in ECI
% outputs a normal vector (1 col) of the magnetic field in ECI frame

% Test r_eci = [7000000;6100000; 5100000]; %m
% t  = 1;
% th_g0 = 175.1465*deg2rad; %greenwhich sidereal time at epoch 

function output = Get_Mag(t,th_g0, r_eci)
deg2rad = pi/180;

r_hat = r_eci/norm(r_eci);  %normal vector of satellite

deg2rad = pi/180;
R_earth = 6378*1000;    %m
H0 = 30115; %nT
th_m = 196.54*deg2rad;  %co-evelation at 2003 (deg)
psi_m = 108.43*deg2rad; %East longitude of dipole at 2003 (deg)
omega = 7.2921159e-5;   %average rotation rate of the earth rad/sec

alpha_m = th_g0 + omega*t + psi_m;

d = [sin(th_m)*cos(alpha_m);
    sin(th_m)*sin(alpha_m);
    cos(th_m)];

%Components of the dipole model in ECI frame
m_i = R_earth^3*H0/(norm(r_eci)^3)*3*d'*r_hat*r_hat-d;

output = m_i/norm(m_i); %normalise output

end