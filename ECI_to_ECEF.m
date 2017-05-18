function [output] = ECI_to_ECEF(input)
% ECI_to_ECEF Converts ECI coordinates to ECEF
%   Take in a vector of the form 
%   input = [x; y; z; t]
%   And outputs a vector in ECEF coordinates of the form 
%   output = [xx; yy; zz]

omega_ie = 7.292115e-5; %rad/sec
%get componenets from input vector
x = input(1);
y = input(2);
z = input(3);
t = input(4);

ECI_to_ECEF_matrix = [
    cos(omega_ie*t),sin(omega_ie*t),0;
    -sin(omega_ie*t),cos(omega_ie*t),0;
    0,  0,  1]; 
output = ECI_to_ECEF_matrix*[x; y; z];

end