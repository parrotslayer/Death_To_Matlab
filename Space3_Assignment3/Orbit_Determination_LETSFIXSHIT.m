%% Assignment 2 Q2
%outputs the ECEF of the satellite and the estimated ECEF of the satellite

clc 
clear all
close all

orbit_params(1)= 7162*1000;              %a - semi major axis meters    
orbit_params(2)=0.0000872;                %e - eccentricity deg
orbit_params(3)=98.7401;                  %inc - inclination degrees
orbit_params(4)=142.1145;                 %Omega - degrees
orbit_params(5)=25.5213;                  %omega - degrees
orbit_params(6)=283.570800000000;         %Mo - Mean Anomaly
orbit_params(7)=2457866.50000000;         %Julian Day (Epoch) Sunday 23/4/17 UT1 00:00:00

%function [Sat_ECI_true,Sat_ECEF_true,Sat_ECI_est, Sat_ECEF_est] = Orbit_Determination(orbit_params)
% Constants
% Define variable u = GM, M is mass of earth
G = 6.67408e-11 ;
massearth = 5.972e24;
u = G * massearth;
rad2deg = 180/pi;
deg2rad = pi/180;

%% Get GLONASS data
PlotEarthLatLong
gs_names = {'Adelaide, South Australia','Hermitage, UK','Cape Caneveral, Florida'};
gs_llh = [-33.9284989  ,  138.6007456   , 45;...
    51.4520259, -1.2754 , 110  ;...
    28.3922182,  -80.6077132 , 2];
% Plot the position of the ground stations
for t = 1 : length(gs_llh(:,1))
    hold on
    plot(gs_llh(t,2),gs_llh(t,1),'MarkerSize',5,'LineWidth',2,'Marker','x','Color','r')
end

% Northern: Monday, 20 March 2017 at 10:29 UTC
julian_date17 = juliandate(datetime([2017,3,20,10,29,0]))*86400;

time_epoch = max(orbit_params(:,7)) * 86400;  %in seconds
time_period = 60 * 60 * 24;  %24 hours

num_sat = 1;  %number of satellites
num_gs = length(gs_llh(:,1));   %number of ground stations

%pre-allocate ECI, ECEF, LGCV, Polar, Time after Epoch arrays filled with NaN
ECI_sat = NaN(3,time_period,num_sat);
ECI_sat_vel = NaN(time_period,1,num_sat);
ECEF_sat = NaN(3,time_period,num_sat);
LGCV_sat = NaN(3,time_period,num_gs);
RAE_sat = NaN(3,time_period,num_gs);
RAE_FULL = NaN(3,time_period,num_gs);
tafterepoch = NaN(num_sat);

%% For each satellite, plot and calculate the ECI data
for k = 1 : num_sat
    % Plot and calculate ECI data starting from the last measured ephemeris
    % Time that has expired until the last measured satellite
%time since Mean anomaly
    tafterepoch(k) = time_epoch - orbit_params(k,7)*86400;
    [ECI_sat(:,:,k),ECI_sat_vel(:,:,k)] = Orbit_to_ECI_and_Simulate(orbit_params(k,1),orbit_params(k,2),orbit_params(k,3),...
        orbit_params(k,4),orbit_params(k,5),orbit_params(k,6),...
        tafterepoch(k),time_period);
end

%% Get the Range Azimuth Elevation over 24 hours
global step
step = 5;  %if we dont want to calculate every second

% Get ECEF coordinates of ALL satellites over 24 hour period
for t = 1:step:time_period
    current_time = time_epoch + t;    %time since last epoch for time = n
    % Time since the last equinox, for ECEF coordinate calculation
    tsinceequinox = (current_time - julian_date17);
    % Get the ECEF coordinates of all satellites at this specific time
    for j = 1 : num_sat
        ECEF_sat(:,t,j) = ECI_to_ECEF([ECI_sat(:,t,j);tsinceequinox]);
    end
end

% Get the LGCV, RAE of a single satellite over 24hours when visible
% Choose satellite number 1
sat_num = 1;
% loop for timesteps in 24h
for t = 1:step:time_period
    %repeat for all ground stations
    for k = 1 : num_gs
        %convert ECEF of Satellite to LGCV
        LGCV_sat(:,t,k) = ECEF_to_LGCV(gs_llh(k,1), gs_llh(k,2), gs_llh(k,3),...
            ECEF_sat(1,t,sat_num), ECEF_sat(2,t,sat_num), ECEF_sat(3,t,sat_num));
        
        %convert LGCV to RAE (degrees)
        RAE_sat(:,t,k) = LGCV_to_RAE_deg(LGCV_sat(:,t,k));
        RAE_FULL(:,t,k) = RAE_sat(:,t,k);
        
        %Set to NaN if Elevation is less than 0 ie below horizon
        if RAE_sat(3,t,k) < 0
            RAE_sat(:,t,k) = NaN;
        end
    end
end

%% Plot the RAE of each satellite when visible and all RAE

figure
subplot(3,1,1)
plot ((1:time_period),RAE_sat(1,:,1),'.b');
hold on
plot ((1:time_period),RAE_sat(1,:,2),'.r');
hold on
plot ((1:time_period),RAE_sat(1,:,3),'.g');
legend('Adelaide, South Australia','Hermitage, UK','Cape Caneveral, Florida' )
title('Elevation vs Time')
xlabel('Time (seconds)')
ylabel('Range (m)')

subplot(3,1,2)
plot ((1:time_period),RAE_sat(2,:,1),'.b');
hold on
plot ((1:time_period),RAE_sat(2,:,2),'.r');
hold on
plot ((1:time_period),RAE_sat(2,:,3),'.g');
title('Azimuth vs Time')
xlabel('Time (seconds)')
ylabel('Azimuth (deg)')

subplot(3,1,3)
plot ((1:time_period),RAE_sat(3,:,1),'.b');
hold on
plot ((1:time_period),RAE_sat(3,:,2),'.r');
hold on
plot ((1:time_period),RAE_sat(3,:,3),'.g');
title('Elevation vs Time')
xlabel('Time (seconds)')
ylabel('Elevation (deg)')
%% Add Tracking Errors to Readings
% LASER RANGING - 1arcsec precision, 5-10mm range error
% Preallocating arrays for speed (TY matlab)
range = zeros(num_gs,time_period);
azimuth = zeros(num_gs,time_period);
elevation = zeros(num_gs,time_period);

gain = 1;
range_error = gain*0.01; %10mm error
angle_error = gain*1/3600;   %1 arc second error

%use gaussian errors
for k = 1:num_gs
range(k,:) = normrnd(RAE_sat(1,:,k), range_error);
azimuth(k,:) = normrnd(RAE_sat(2,:,k), angle_error);
elevation(k,:) = normrnd(RAE_sat(3,:,k), angle_error);
end

%% Initial Estimate of R_ECI and V_ECI
% From the average of guesses from first 3 observations of each,
% Utilising herricks gibbs to get estimated velocities

%Pre-allocate arrays of NaN to be filled with real measurements
init_guess_input = NaN(4,3,num_gs);
indexes = NaN(num_gs,3);
store_init_guess = NaN(6,num_gs);
nlls_orbit = NaN(5,6,3);
params_error = NaN(5,6,3);

for k = 1 : num_gs
    % get ALL the indexes of RAE entries that are not NaN i.e. visible
    indexes = find(not(isnan(azimuth(k,:))));
    
    increment = 200;
    for j = 1:length(indexes)/increment
    start_i = (j-1)*increment+1;
    end_i = j*increment+1;
    
        % Take the first 3 observations from each station with real values
    indexes_obs(k,:) = indexes(start_i:start_i+3);
    
    %get time of second measurement for each ground station
    time2(j,k) = indexes_obs(k,2);
    
    %Uses Herricks Gibbs to estimate the velocity using ECI measurements
    %and to also return the ECI position corresponding to the velocity
    %estimate (2nd measurement)
    init_guess(:,k) = get_init_posvel([time_epoch+indexes_obs(k,:);
        range(k,indexes_obs(k,:));
        azimuth(k,indexes_obs(k,:));
        elevation(k,indexes_obs(k,:))],...
        julian_date17,...
        gs_llh(k,:));...
    
    % Take all observations when the satellite is visible
    obs = [time_epoch + indexes(start_i:end_i);...
        range(k,indexes(start_i:end_i));...
        azimuth(k,indexes(start_i:end_i))*pi/180;...
        elevation(k,indexes(start_i:end_i))*pi/180]';
    
    posvel_guess = init_guess(:,k);
    
    %convert ground staiton LLH to ECEF
    GS_ECEF(1,:) = LLH_to_ECEF_Geocentric(gs_llh(k,1), gs_llh(k,2), gs_llh(k,3));
    
    %run nlls function given from lecturer
    % output = Orbital parameters: [a,e,i,RAAN,AoP,Mo(at epoch)] m, degrees
    nlls_orbit(j,:,k) = nlls_orbit_determ(julian_date17,obs,GS_ECEF,posvel_guess,gs_llh(k,:));
    
    error_a = nlls_orbit(j,1,k) - orbit_params(1);   %in m
    error_e = nlls_orbit(j,2,k) - orbit_params(2);
    error_i = nlls_orbit(j,3,k) - orbit_params(3);
    error_Omega = nlls_orbit(j,4,k) - orbit_params(4);
    error_omega = nlls_orbit(j,5,k) - orbit_params(5);
    
    mu = 3.986005e14;
    n = sqrt(mu/(orbit_params(1,1))^3);
    MeanAnomaly = orbit_params(1,6) + n*(time2(j,k)*deg2rad);    %degrees
    error_M = nlls_orbit(j,6,k) - MeanAnomaly;  %mean anomaly is a function of time
    
    params_error(j,:,k) = [error_a, error_e, error_i, error_Omega, error_omega, error_M];
    
    end
end


%% Get Estimated ECI Orbit using orbital parameters from NLLS for a useful ground station
%some will fail to converge and give NaN so use one that is good.
global gs_num
% for k = 1:num_gs
%     if not(isnan(nlls_orbit(k,1)))
%         gs_num = k;
%         disp(['Using Ground Station: ' gs_names{k}]);
%         break;
%     end
% end
gs_num = 2;
j = 1;
a_est = nlls_orbit(j,1,gs_num);
e_est = nlls_orbit(j,2,gs_num);
i_est = nlls_orbit(j,3,gs_num);
Omega_est = nlls_orbit(j,4,gs_num);
omega_est = nlls_orbit(j,5,gs_num);
M_est = nlls_orbit(j,6,gs_num);

ff = 227;   %offset to sync estimated and real ECI values 
%ff = -401;
ff = 0;

% Need mean anomaly at the time at which herricks gibbs was taken
mu = 3.986005e14;
n = sqrt(mu/(orbit_params(1,1))^3);
MeanAnomaly = orbit_params(1,6) + n*(time2(j,gs_num))*deg2rad;
k = 1;
[ECI_sat2(:,:),ECI_sat_vel2(:,:)] = Orbit_to_ECI_and_Simulate(orbit_params(k,1)...
    ,orbit_params(k,2),orbit_params(k,3),orbit_params(k,4),orbit_params(k,5),...
    MeanAnomaly,0,time_period);

%get estimated eci orbit from NLLS estimation of a single ground station
ECI_est = Orbit_to_ECI_and_Simulate(a_est,e_est,i_est,Omega_est,omega_est,...
    MeanAnomaly,0,time_period);
%% error plots

% get error for the specified satellite
% error_a = a_est - orbit_params(1);   %in m
% error_e = e_est - orbit_params(2);
% error_i = i_est - orbit_params(3);
% error_Omega = Omega_est - orbit_params(4);
% error_omega = omega_est - orbit_params(5);
% error_M = M_est - orbit_params(6);
% 
% params_error = [error_a, error_e, error_i, error_Omega, error_omega, error_M];

error = NaN(time_period,1);

for t = 1:time_period
    dx = ECI_est(1,t) - ECI_sat2(1,t);
    dy = ECI_est(2,t) - ECI_sat2(2,t);
    dz = ECI_est(3,t) - ECI_sat2(3,t);
    error(t) = sqrt(dx^2 + dy^2 + dz^2);
end

figure
plot(1:time_period,error)
xlabel('Time (seconds)')
ylabel('Position Error (m)')
title('Position Error vs Time')

%% Plot ECI estimated and Real X position vs time
figure
plot(1:time_period,ECI_est(1,:),'r')
hold on
plot(1:time_period,ECI_sat2(1,:),'b')
legend('Estimated X Position','Real X Position')
title('ECI X position vs Time')
xlabel('Time (seconds)')
ylabel('Position (m)')

%% convert estimated ECI orbit to ECEF
% Get ECEF and LGCV coordinates from the ECI coordinates
for t = 1:step:time_period
    current_time = time_epoch + t;    %time since last epoch for time = n
    % Time since the last equinox, for ECEF coordinate calculation
    tsinceequinox = (current_time - julian_date17);
    % Get the ECEF coordinates of all satellites at this specific time
    ECEF_est(:,t) = ECI_to_ECEF([ECI_sat(:,t,j);tsinceequinox]);
    % Get ECEF coordinates wrt the chosen ground station
    LGCV_est(:,t) = ECEF_to_LGCV(gs_llh(gs_num,1), gs_llh(gs_num,2), gs_llh(gs_num,3),...
        ECEF_sat(1,t), ECEF_sat(2,t), ECEF_sat(3,t));
end

%% Plot ECI estimated and real over 24h
% plot the 3D ECI Orbits from the 1st ground measurement

PlotEarthSphere
hold on
plot3(ECI_sat(1,:,sat_num),ECI_sat(2,:,sat_num),ECI_sat(3,:,sat_num),'b')
hold on
plot3(ECI_est(1,:),ECI_est(2,:),ECI_est(3,:),'r.')
legend('Earth','Real Orbit','Estimated Orbit')

%% Output variables
%Sat_ECI_true,Sat_ECEF_true,Sat_ECI_est, Sat_ECEF_est
Sat_ECI_true = ECI_sat;
Sat_ECEF_true = ECEF_sat(:,:,gs_num);
Sat_ECI_est = ECI_est;
Sat_ECEF_est = ECEF_est;

%end