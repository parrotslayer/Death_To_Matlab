% Calculates required FOV for a given Swath width and Height

H = 7.8335e5;   %m
Swath_Width = 200e3;  %m

n = atan((Swath_Width/2)/H);

FOV_rads = 2*n  %radians

FOV_deg = FOV_rads*180/pi