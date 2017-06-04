% Converts orbital parameters to ECI and then simulates it over a
% timeperiod
%input units in m and degrees and seconds
function [P_ECI,velocity] = Orbit_to_ECI_and_Simulate(a,e,i,omega,w,Mo,tsinceMo,timeperiod)
% a,e,i,w,omega,Mo,time since measured parameters, period of time to
% simulate over
deg2rad = pi/180 ;
i = deg2rad*i;
w = deg2rad*w;
omega = deg2rad*omega;
Mo = deg2rad*Mo;
mu = 3.986005e14;

% Transformation matrix C orbit -> ECI
C(1,:) = [cos(omega)*cos(w) - sin(omega)*sin(w)*cos(i),...
    -cos(omega)*sin(w) - sin(omega)*cos(w)*cos(i),...
    sin(omega)*sin(i)];
C(2,:) = [sin(omega)*cos(w) + cos(omega)*sin(w)*cos(i),...
    -sin(omega)*sin(w) + cos(omega)*cos(w)*cos(i),...
    -cos(omega)*sin(i)];
C(3,:) = [sin(w)*sin(i),...
    cos(w)*sin(i),...
    cos(i)];

% Mean motion of orbit
n = sqrt(mu/a^3);

% Semi-lactus rectum of orbit
p = a*(1-e^2);

%Pre-allocate for speed
E_matrix = zeros(timeperiod,1);
theta = zeros(timeperiod,1);
r = zeros(timeperiod,1);
P_ECI = zeros(3,timeperiod,1);
velocity = zeros(timeperiod,1);

%% Simulate ECI over a timeperiod
for k = 1:timeperiod
    % Solving for E , M = E - e*sin(E)
    %M = wrapTo2Pi(Mo + (k+tsinceMo)*n);
    M = Mo + (k+tsinceMo)*n;
    E = M;
    
    % Define the function for Newton-Raphson method
    func = E - e*sin(E) - M;
    % Iterate until convergence condition is met using Newton-Raphson
    thresh = 0.001;     %threshold for convergence
    while abs(func) > thresh
        E = E - (E - e*sin(E)- M)/(1-e*cos(E));
        func = E - e*sin(E) - M;
    end
    
    % Store value of Eccentric anomaly
    E_matrix(k) = E;
    
    % Solving Kepler's Equation for true anomaly from eccentric
    theta(k) = 2*atan2(((1+e)/(1-e))^(0.5)*sin(E/2) , cos(E/2));
    
    % Solve for radius given theta
    r(k) = p/(1+e*cos(theta(k)));
    
    % Orbital coordinates
    orbit_x = r(k)*cos(theta(k)); %Add intitial starting angle
    orbit_y = r(k)*sin(theta(k)); 
    orbit_z = 0 ;
    orbit_vec = [orbit_x;orbit_y;orbit_z];
    %apply conversion matrix
    P_ECI(:,k) = C*orbit_vec;
    
    % Calculate Velocity at the Time
    velocity(k) = sqrt(mu*(2/r(k) - 1/a));
end
end