% Build C_Body_to_LGCV
% Takes in [Yaw, Pitch, Roll] in radians
% [Psi, Th, Phi]
function C = C_LGCV_to_Body(attitude)
psi = attitude(1);
th = attitude(2);  
Phi = attitude(3);  

C(1,:) = [cos(psi)*cos(th), sin(psi)*cos(th), -sin(th)];
C(2,:) = [(cos(psi)*sin(th)*sin(Phi)-sin(psi)*cos(Phi)),...
    (sin(psi)*sin(th)*sin(Phi)+cos(psi)*cos(Phi)), cos(th)*sin(Phi) ];
C(3,:) = [(cos(psi)*sin(th)*cos(Phi)+sin(psi)*sin(Phi)),...
    (sin(psi)*sin(th)*cos(Phi)-cos(psi)*sin(Phi)), cos(th)*cos(Phi)];
end