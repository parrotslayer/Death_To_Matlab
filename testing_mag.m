%% Get Magnetometer and Star Tracker Readings Over 24h
clc
th_g0 = 175.1465*deg2rad; %greenwhich sidereal time at epoch
%get the LLH of the chosen GS
chosen_gs_LLH = gs_llh(chosen_gs,:);
time_epoch = orbit_params(7) * 86400;  %time of epoch in seconds


for t = 1:step:time_period
 % Test r_eci = [7000000;6100000; 5100000]; %m
Mag_ECI(t,:) = Get_Mag(t,th_g0, ECI_sat_est(:,t));

% Convert ECI to ECEF
current_time = time_epoch + t;    %time since last epoch for time = n
t_since_equinox = (current_time - julian_date17);
Mag_ECEF(t,:) = ECI_to_ECEF([ECI_sat_est(:,t); t_since_equinox]);

Mag_LGCV(t,:) =  ECEF_to_LGCV(chosen_gs_LLH(1),chosen_gs_LLH(2),chosen_gs_LLH(3),...
    Mag_ECEF(t,1), Mag_ECEF(t,1), Mag_ECEF(t,1));  
end