% lambda_i and phi are in degrees
% R in meters

function [output] = LLH_to_ECEF_Geodetic(lambda, phi, height)
radius_earth = 6378137; %from lecture slides
e = 0.081819;
N = radius_earth/( sqrt(1-e^2*sind(lambda)^2) ) ;
rx = (N+height)*cosd(lambda)*cosd(phi);
ry = (N+height)*cosd(lambda)*sind(phi); 
rz = (N*(1-e^2)+height)*sind(lambda); 
output = [rx, ry, rz];
end