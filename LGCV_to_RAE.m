%% LGCV to RAE 
% local geocentric vertical to range azimuth (phi) elevation (theta)
% input = [x; y; z] in meters
% output = [range; azimuth, elevation] in m and radians
function [output] = LGCV_to_RAE(input)

x = input(1);
y = input(2);
z = input(3);

R = sqrt(x^2 + y^2 + z^2);
phi = atan2(y,x);
theta = atan( -z/sqrt(x^2+y^2) );

output = [R;phi;theta];

end