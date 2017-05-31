% Rotates a 3D vector by the 3 euler angles 
% Takes in [Roll, Pitch, Yaw] [Phi, Th, psi] in radians
% Takes in a column vector [X,Y,Z]
% Returns a column vector that is rotated XYZ
% CONFIRMED

function output = Rotate_3D(angles,input)
Phi = angles(1);  %actually roll
th = angles(2);   %pitch
psi = angles(3);  %actually yaw

C(1,:) = [cos(psi)*cos(th), sin(psi)*cos(th), -sin(th)];
C(2,:) = [(cos(psi)*sin(th)*sin(Phi)-sin(psi)*cos(Phi)),...
    (sin(psi)*sin(th)*sin(Phi)+cos(psi)*cos(Phi)), cos(th)*sin(Phi) ];
C(3,:) = [(cos(psi)*sin(th)*cos(Phi)+sin(psi)*sin(Phi)),...
    (sin(psi)*sin(th)*cos(Phi)-cos(psi)*sin(Phi)), cos(th)*cos(Phi)];

%return a column vector
output = (C*input')';

end