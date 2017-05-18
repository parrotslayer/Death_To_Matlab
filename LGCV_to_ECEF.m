%% LGCV to ECEF
%converts LGCV coordinates + LLH of ground station to ECEF
%input = lat_g, long_g, height_g, rx, ry, rz
%output = [x; y; z]

function R_ECEF = LGCV_to_ECEF(lat_g, long_g, height_g, rx, ry, rz)
%conversion matrix
ECEF_to_LGCV_matrix = [
    -sind(lat_g)*cosd(long_g),  -sind(long_g),  -cosd(lat_g)*cosd(long_g);
    -sind(lat_g)*sind(long_g),  cosd(long_g),   -cosd(lat_g)*sind(long_g);
    cosd(lat_g), 0,  -sind(lat_g)];

% R_ECEF = [matrix]*P_LGCV
X_rel = ECEF_to_LGCV_matrix*[rx; ry; rz];

%get observer location in ECEF
X_obs = transpose(LLH_to_ECEF(lat_g,long_g,height_g));

%get absolute vector from the origin
R_ECEF = X_rel + X_obs;
end