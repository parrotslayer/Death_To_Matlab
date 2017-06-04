% Build Stars
% Takes in orbital parameters [inc, Omega, omega] in radians
% as well as r = radius of stars and the number of stars to plot
% Returns star's ECI coordinates
function output = Summon_Stars(r,inc,Omega,omega,num_stars)
theta = linspace(0,2*pi,num_stars);    %creates 1000 equally spaced components from 0 to 2pi

% Transformation matrix C orbit -> ECI
C(1,:) = [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(inc),...
    -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(inc),...
    sin(Omega)*sin(inc)];
C(2,:) = [sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(inc),...
    -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(inc),...
    -cos(Omega)*sin(inc)];
C(3,:) = [sin(omega)*sin(inc),...
    cos(omega)*sin(inc),...
    cos(inc)];

% Orbital coordinates to be outputted
for i = 1:num_stars
    orbit_x = r*cos(theta(i));
    orbit_y = r*sin(theta(i));
    orbit_z = 0;
    orbit = [orbit_x;orbit_y;orbit_z];
    Star_ECI(i,:) = C*orbit;
end
output = Star_ECI;

end
