%% Polar to Cart  
% range azimuth (phi) elevation (theta) to cartesian 
% input = [range; azimuth, elevation] in m and rads
% output = [x; y; z] % in m
function [output] = Polar_to_Cart(input)

R = input(1);
phi = input(2);
theta = input(3);

%in degrees so use cosd
x = R*cos(theta)*cos(phi);
y = R*cos(theta)*sin(phi);
z = -R*sin(theta);

output = [x;y;z];

end