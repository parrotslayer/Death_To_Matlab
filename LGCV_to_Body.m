% Build Body_to_LGCV
% Takes in [Yaw, Pitch, Roll] [Phi, Th, psi] in radians
% Takes in LGCV coordinates as a column vector [X,Y,Z]
% Returns a column vector
% CONFIRMED

function output = LGCV_to_Body(attitude,LGCV)
Phi = attitude(1);  %yaw
th = attitude(2);   %pitch
psi = attitude(3);  %roll

C(1,:) = [cos(psi)*cos(th), sin(psi)*cos(th), -sin(th)];
C(2,:) = [(cos(psi)*sin(th)*sin(Phi)-sin(psi)*cos(Phi)),...
    (sin(psi)*sin(th)*sin(Phi)+cos(psi)*cos(Phi)), cos(th)*sin(Phi) ];
C(3,:) = [(cos(psi)*sin(th)*cos(Phi)+sin(psi)*sin(Phi)),...
    (sin(psi)*sin(th)*cos(Phi)-cos(psi)*sin(Phi)), cos(th)*cos(Phi)];

%return a column vector
output = (C*LGCV')';

end