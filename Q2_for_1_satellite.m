%% Assignment 2 Q2
clear
clc
close all

function [output1,output2] = Get_Estimated_ECEF()

% Constants
% Define variable u = GM, M is mass of earth
G = 6.67408e-11 ;
massearth = 5.972e24;
u = G * massearth;
rad2deg = 180/pi;
deg2rad = pi/180;

%% Get GLONASS data

PlotEarthLatLong
% ground station names and locations
gs_names = {'Schyolkovo','Komsomolsk','St.Petersburg'};
gs_llh = [55.910087  ,  38.009208   , 150;...
    49.0156679, 33.6450451 , 72  ;...
    59.9342802,  30.3350986 , 17];

% Plot the position of the ground stations
for t = 1 : length(gs_llh(:,1))
    hold on
    plot(gs_llh(t,2),gs_llh(t,1),'MarkerSize',5,'LineWidth',2,'Marker','x','Color','r')
end

%get the GLONASS data from a MAT file which converted TLE data already
load GLONASS.mat
%satellite ephermis data
ephsat = GLONASS(1,(2:end));

ephsat(1) = 7162*1000;              %a - semi major axis meters    
ephsat(2)=0.0000872;                %e - eccentricity deg
ephsat(3)=98.7401;                  %inc - inclination degrees
ephsat(4)=142.1145;                 %Omega - degrees
ephsat(5)=25.5213;                  %omega - degrees
ephsat(6)=283.570800000000;         %Mo - Mean Anomaly
ephsat(7)=2457866.50000000;         %Time since epoch

% Northern: Monday, 20 March 2017 at 10:29 UTC
julian_date17 = juliandate(datetime([2017,3,20,10,29,0]))*86400;

time_epoch = max(ephsat(:,7)) * 86400;  %in seconds
time_period = 60 * 60 * 24;  %24 hours

num_sat = 1;  %number of satellites
num_gs = length(gs_llh(:,1));   %number of ground stations

%pre-allocate ECI, ECEF, LGCV, Polar, Time after Epoch arrays filled with NaN
ECI_sat = NaN(3,time_period,num_sat);
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
    tafterepoch(k) = time_epoch - ephsat(k,7)*86400;
    ECI_sat(:,:,k) = Orbit_to_ECI_and_Simulate(ephsat(k,1),ephsat(k,2),ephsat(k,3),...
        ephsat(k,4),ephsat(k,5),ephsat(k,6),...
        tafterepoch(k),time_period);
end

% plot the 3D ECI Orbits
PlotEarthSphere
hold on
for k = 1 : num_sat
    plot3(ECI_sat(1,:,k),ECI_sat(2,:,k),ECI_sat(3,:,k))
end

%% Get the Range Azimuth Elevation over 24 hours
global step
step = 100;  %if we dont want to calculate every second

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
% issue with step because if step ~= 1 the plot for all RAE wont show
if step == 100
    for k = 1 : num_gs
        string = gs_names{k};
        figure
        subplot(3,1,1)
        plot ((1:time_period),RAE_sat(1,:,k),'.b');
        hold on
        plot ((1:time_period),RAE_FULL(1,:,k),'r');
        title ({string;'Range (m)'})
        
        subplot(3,1,2)
        plot ((1:time_period),RAE_sat(2,:,k),'.b');
        hold on
        plot ((1:time_period),RAE_FULL(2,:,k),'r');
        title ('Azimuth (deg)')
        
        subplot(3,1,3)
        plot ((1:time_period),RAE_sat(3,:,k),'.b');
        hold on
        plot ((1:time_period),RAE_FULL(3,:,k),'r');
        title ( 'Theta (Elevation) (deg) ' )
        xlabel('Time (s)')
        legend('When Visible', 'Entire Duration')
    end
else
    for k = 1 : num_gs
        string = gs_names{k};
        figure
        subplot(3,1,1)
        plot ((1:time_period),RAE_sat(1,:,k),'.b');
        title ({string;'Range (m)'})
        
        subplot(3,1,2)
        plot ((1:time_period),RAE_sat(2,:,k),'.b');
        title ('Azimuth (deg)')
        
        subplot(3,1,3)
        plot ((1:time_period),RAE_sat(3,:,k),'.b');
        title ( 'Theta (Elevation) (deg) ' )
        xlabel('Time (s)')
    end
end

%% Add Tracking Errors to Readings
% LASER RANGING - 1arcsec precision, 5-10mm range error
% Preallocating arrays for speed (TY matlab)
range = zeros(num_gs,time_period);
azimuth = zeros(num_gs,time_period);
elevation = zeros(num_gs,time_period);

gain = 1;
range_error = gain*0.01; %10mm error
angle_error = gain*1/3600;   %1 arc second error

for t = 1 : time_period
    for k = 1 : num_gs
%         range(k,t) = RAE_sat(1,t,k); + range_error*(-1 + rand);
%         azimuth(k,t) = RAE_sat(2,t,k); + angle_error*(-1 + rand);
%         elevation(k,t) = RAE_sat(3,t,k) + angle_error*(-1 + rand);

%use gaussian errors
        range(k,t) = normrnd(RAE_sat(1,t,k), range_error);
        azimuth(k,t) = normrnd(RAE_sat(2,t,k), angle_error);
        elevation(k,t) = normrnd(RAE_sat(3,t,k), angle_error);
    end
end

%% Initial Estimate of R_ECI and V_ECI
% From the average of guesses from first 3 observations of each,
% Utilising herricks gibbs to get estimated velocities

%Pre-allocate arrays of NaN to be filled with real measurements
init_guess_input = NaN(4,3,num_gs);
indexes = NaN(num_gs,3);
store_init_guess = NaN(6,num_gs);

for k = 1 : num_gs
    % Take the first 3 observations from each station with real values
    indexes(k,:) = find(not(isnan(azimuth(k,:))),3);
    
    %Uses Herricks Gibbs to estimate the velocity using ECI measurements
    %and to also return the ECI position corresponding to the velocity
    %estimate (2nd measurement)
    init_guess(:,k) = get_init_posvel([time_epoch+indexes(k,:);
        range(k,indexes(k,:));
        azimuth(k,indexes(k,:));
        elevation(k,indexes(k,:))],...
        julian_date17,...
        gs_llh(k,:));...
end

%% Implement NLLS Universal Conic Section on Satellite 1
for k = 1 : num_gs
    % get the indexes of RAE entries that are not NaN i.e. visible
    indexes = find(not(isnan(azimuth(k,:))));
    
    % Take all observations when the satellite is visible
    obs = [time_epoch + indexes(1:end);...
        range(k,indexes(1:end));...
        azimuth(k,indexes(1:end))*pi/180;...
        elevation(k,indexes(1:end))*pi/180]';
    
    posvel_guess = init_guess(:,k);
    
    %convert ground staiton LLH to ECEF
    %GS_ECEF(1,:) = gcLLHtoECEF(gs_llh(k,:));
    GS_ECEF(1,:) = LLH_to_ECEF_Geocentric(gs_llh(k,1), gs_llh(k,2), gs_llh(k,3));
    
    %run nlls function given from lecturer
    % output = Orbital parameters: [a,e,i,RAAN,AoP,Mo(at epoch)] m, degrees
    nlls_orbit(k,:) = nlls_orbit_determ(julian_date17,obs,GS_ECEF,posvel_guess,gs_llh(k,:));
end

%% Get Estimated ECI Orbit using orbital parameters from NLLS for a useful ground station
%some will fail to converge and give NaN so use one that is good.
for k = 1:num_gs
    if not(isnan(nlls_orbit(k,1)));
        gs_num = k;
        disp(['Using Ground Station: ' gs_names{k}]);
        break;
    end
end
 
a_est = nlls_orbit(gs_num,1);
e_est = nlls_orbit(gs_num,2);
i_est = nlls_orbit(gs_num,3);
Omega_est = nlls_orbit(gs_num,4);
omega_est = nlls_orbit(gs_num,5);
M_est = nlls_orbit(gs_num,6);

%time of observation from herrick's gibbs
global time2

%get estimated eci orbit
ECI_est = Orbit_to_ECI_and_Simulate(a_est,e_est,i_est,Omega_est,omega_est,...
    M_est,time2,time_period);

%convert estimated ECI orbit to ECEF
% Get ECEF coordinates of ALL satellites over 24 hour period
for t = 1:step:time_period
    current_time = time_epoch + t;    %time since last epoch for time = n
    % Time since the last equinox, for ECEF coordinate calculation
    tsinceequinox = (current_time - julian_date17);
    % Get the ECEF coordinates of all satellites at this specific time
    ECEF_est(:,t) = ECI_to_ECEF([ECI_sat(:,t,j);tsinceequinox]);
end

output1 = ECEF_sat;
output2 = ECEF_est;
%% Plot ECI estimated and real over 24h
% plot the 3D ECI Orbits from the 1st ground measurement

PlotEarthSphere
hold on
plot3(ECI_sat(1,:,sat_num),ECI_sat(2,:,sat_num),ECI_sat(3,:,sat_num),'b')
hold on
plot3(ECI_est(1,:),ECI_est(2,:),ECI_est(3,:),'r.')
legend('Earth','Real Orbit','Estimated Orbit')