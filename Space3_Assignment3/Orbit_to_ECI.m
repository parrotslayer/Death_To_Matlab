function output = Orbit_to_ECI(a,e,i,Omega,omega,Mo,tsinceMo)
% a,e,i,Omega,omega,Mo,time since measured parameters, period of time to
% units are in m and degrees and seconds

%convert units to radians
deg2rad = pi/180 ;
i = deg2rad*i;
omega = deg2rad*omega;
Omega = deg2rad*Omega;
Mo = deg2rad*Mo;

% Transformation matrix C orbit -> ECI
C(1,:) = [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(i),...
        -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(i),...
        sin(Omega)*sin(i)];
C(2,:) = [sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(i),...
        -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(i),...
        -cos(Omega)*sin(i)];
C(3,:) = [sin(omega)*sin(i),...
         cos(omega)*sin(i),...
         cos(i)];

mu = 3.986005e14;    %mu from assignment

% Mean motion of orbit
n = sqrt(mu/a^3);

% Semi-lactus rectum (p) of orbit
p = a* (1-e^2);

% Solving for E , M = E - e*sin(E)
M = wrapTo2Pi(Mo + tsinceMo*n);
E = M;

% Define the function for Newton-Raphson method
func = E - e*sin(E) - M;

% Iterate until convergence condition is met using Newton-Raphson
tol_E = 0.00001;    %tolerance for convergence

while abs(func) > tol_E
E = E - (E - e*sin(E)- M)/(1-e*cos(E));
func = E - e*sin(E) - M;
end

% Solving Kepler's Equation for true anomaly from eccentric
theta = 2*atan2(((1+e)/(1-e))^(0.5)*sin(E/2) , cos(E/2));

% Solve for radius given theta
r = p/(1+e*cos(theta));

% Orbital coordinates to be outputted
orbit_x = r*cos(theta); 
orbit_y = r*sin(theta); 
orbit_z = 0 ;
orbit = [orbit_x;orbit_y;orbit_z];
ECI_xfer = C*orbit;
output = ECI_xfer;
end