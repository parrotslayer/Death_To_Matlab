% Takes in the attitude [Yaw, Pitch, Roll] and the 3D LGCV measurement
% and outputs a 3x3 jacobian for this measurement

%CONFIRMED in NLLS_attitude_WORKING
function H = Jacobian_Attitude(attitude, LGCV)
Phi = attitude(1);
th = attitude(2);
psi = attitude(3);
mx = LGCV(1);
my = LGCV(2);
mz = LGCV(3);

H(1,1) = 0;
H(1,2) = -cos(psi)*sin(th)*mx - sin(psi)*sin(th)*my - cos(th)*mz;
H(1,3) = -sin(psi)*cos(th)*mx + cos(psi)*cos(th)*my;

H(2,1) = (cos(psi)*sin(th)*cos(Phi)+sin(psi)*sin(Phi))*mx...
    + (sin(psi)*sin(th)*sin(Phi)-cos(psi)*sin(Phi))*my + cos(th)*cos(Phi)*mz;
H(2,2) = cos(psi)*cos(th)*sin(Phi)*mx + sin(Phi)*cos(th)*sin(Phi)*my...
    - sin(th)*sin(Phi)*mz;
H(2,3) = (-sin(psi)*sin(th)*sin(Phi)-cos(psi)*cos(Phi))*mx...
    + (cos(psi)*sin(th)*sin(Phi)-sin(psi)*cos(th))*my;

H(3,1) = (-cos(psi)*sin(th)*sin(Phi) +sin(psi)*cos(Phi))*mx...
    + (-sin(psi)*sin(th)*sin(Phi)-cos(psi)*cos(Phi))*my - cos(th)*sin(Phi)*mz;
H(3,2) = cos(psi)*cos(th)*cos(Phi)*mx + sin(psi)*cos(th)*cos(Phi)*my...
    - sin(th)*cos(Phi)*mz;
H(3,3) = (-sin(psi)*sin(th)*cos(Phi)+cos(psi)*sin(Phi))*mx...
    + (cos(psi)*sin(th)*cos(Phi)+sin(psi)*sin(Phi))*my;

end