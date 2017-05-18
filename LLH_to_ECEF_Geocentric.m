% lambda_i and phi are in degrees
% R in meters

function [output] = LLH_to_ECEF_Geocentric(lambda_i, phi, height)
radius_earth = 6378137; %from lecture slides
R = radius_earth + height;
rx = R*cosd(lambda_i)*cosd(phi);
ry = R*cosd(lambda_i)*sind(phi); 
rz = R*sind(lambda_i); 
output = [rx, ry, rz];
end