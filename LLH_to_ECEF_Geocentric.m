% lambda_i and phi are in degrees
% height above ground level in meters

% CONFIRMED
% Test LLH_to_ECEF_Geocentric(-33, 151, 52)
% Output = -4678515; 2593343, -3473810
function [output] = LLH_to_ECEF_Geocentric(lambda_i, phi, height)
radius_earth = 6378137; %from lecture slides
R = radius_earth + height;
rx = R*cosd(lambda_i)*cosd(phi);
ry = R*cosd(lambda_i)*sind(phi); 
rz = R*sind(lambda_i); 
output = [rx, ry, rz];
end