function output = ECEF_to_LGCV(lat_g, long_g, height_g, rx, ry, rz)
% convert LLH ground station to ECEF geocentric of ground station
% units in degrees and meters
% Test = ECEF_to_LGCV(-33,151,52,-4678800,253900,-3474600)

ECEF_gnd = LLH_to_ECEF_Geocentric(lat_g, long_g ,height_g);

%find the relative ECEF between the object and the ground station
ECEF_sat = [rx, ry, rz];
ECEF_rel = ECEF_sat - ECEF_gnd;
x = ECEF_rel(1);
y = ECEF_rel(2);
z = ECEF_rel(3);

ECEF_to_LGCV_matrix = [
    -sind(lat_g)*cosd(long_g),  -sind(long_g),  -cosd(lat_g)*cosd(long_g);
    -sind(lat_g)*sind(long_g),  cosd(long_g),   -cosd(lat_g)*sind(long_g);
    cosd(lat_g), 0,  -sind(lat_g)];
% LGCV = [matrix]^-1 * ECEF
output = ECEF_to_LGCV_matrix\[x; y; z];
end