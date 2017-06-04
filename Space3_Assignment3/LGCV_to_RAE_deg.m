%% LGCV to RAE deg
% local geocentric vertical to range azimuth (phi) elevation (theta)
% input = [x; y; z] in meters
% output = [range; azimuth, elevation] in m and degrees
function output_vec = cartesian2polar(input_vec)
    %test
    %answer = [1000,-140,60];
    %input_vec = [-383.02;-321.39;-866.03];
    x = input_vec(1);
    y = input_vec(2);
    z = input_vec(3);
    R = norm(input_vec);
    azimuth = atan2d(y,x);
    %elevation
    theta = atand(-z/sqrt(x^2+y^2));
    output_vec = [R ; azimuth ; theta];
end