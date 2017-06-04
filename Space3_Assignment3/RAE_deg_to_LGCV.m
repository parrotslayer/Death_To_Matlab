%% RAE to LGCV deg 
% range azimuth (phi) elevation (theta) to local geocentric vertical 
% input = [range; azimuth, elevation] in m and degrees
% output = [x; y; z] % in m
function output = RAE_deg_to_LGCV(input)
%test
%input = [1000,-140,60];
%answer = [-383.02;-321.39;-866.03]
deg2rad = pi/180;
R = input(1);
azimuth = input(2)*deg2rad;
theta = input(3)*deg2rad;
x = R * cos(theta) * cos(azimuth);
y = R * cos(theta) * sin(azimuth);
z = -R * sin(theta);
output = [x;y;z];
end