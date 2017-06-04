function [output] = ECEF_to_ECI(input)
% Converts ECEF coordinates to ECI
%   Take in a vector of the form 
%   input = [x; y; z; t]
%   And outputs a vector in ECI coordinates of the form 
%   output = [xx; yy; zz]
%   CONFIRMED

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
output = inv(ECI_to_ECEF_matrix)*[x; y; z];

end