%% Get Magnetometer and Star Tracker Readings Over 24h
clc
th_g0 = 175.1465*deg2rad; %greenwhich sidereal time at epoch
time_epoch = orbit_params(7) * 86400;  %time of epoch in seconds

%constants for generating constellation
r = 10e10;      %some big number 
inc = orbit_params(3);
Omega = orbit_params(4);
omega = orbit_params(5);
num_stars = 100;
% Generate constellation of stars in ECI
Star_ECI = Summon_Stars(r,inc,Omega,omega,num_stars);
%FOV
FOV = 90*deg2rad;

for t = 1:step:time_period
% Convert satellite's estimated ECEF position to est LLH
Sat_LLH(t,:) = ECEF_to_LLH(Sat_ECEF_est(1,t),Sat_ECEF_est(2,t),Sat_ECEF_est(3,t));

% Get unit vector of magnetometer readings
Mag_ECI(:,t) = Get_Mag(t,th_g0, Sat_ECI_est(:,t));
% Convert ECI to ECEF
current_time = time_epoch + t;    %time since last epoch for time = n
t_since_equinox = current_time - julian_date17;
Mag_ECEF(t,:) = ECI_to_ECEF([Mag_ECI(:,t); t_since_equinox]);

% Convert ECEF to LGCV wrt Satellite's position
Mag_LGCV(t,:) =  ECEF_to_LGCV(Sat_LLH(1),Sat_LLH(2),Sat_LLH(3),...
    Mag_ECEF(t,1), Mag_ECEF(t,1), Mag_ECEF(t,1));  

visible_stars = Get_Star_Tracker(Sat_LLH,t_since_equinox,Star_ECI,FOV);

end