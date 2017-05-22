% NLLS of attitude
clc
close all
clear all

time_period = 24*60*60; %1day
global step
step = 100;         %time step for looping calculations
deg2rad = pi/180;
rad2deg = 180/pi;
% Northern: Monday, 20 March 2017 at 10:29 UTC
julian_date17 = juliandate(datetime([2017,3,20,10,29,0]))*86400;

%% Get Estimated and Real Satellite positions over 24h
orbit_params(1)= 7162*1000;              %a - semi major axis meters    
orbit_params(2)=0.0000872;                %e - eccentricity deg
orbit_params(3)=98.7401;                  %inc - inclination degrees
orbit_params(4)=142.1145;                 %Omega - degrees
orbit_params(5)=25.5213;                  %omega - degrees
orbit_params(6)=283.570800000000;         %Mo - Mean Anomaly
orbit_params(7)=2457866.50000000;         %Julian Day (Epoch)

% Use code from Assignment 2 to model estimated orbit and real orbit
%[Sat_ECEF_true,Sat_ECI_est, Sat_ECEF_est,Sat_LGCV_est] = Get_Estimated_ECEF(orbit_params);
load('Get_Estimated_ECEF_data')

%% Simulate the Attitude
t = 1:time_period;   %generate times
omega1 = 0.001;     %freq of pitch
omega2 = 0.005;     %freq of yaw

%model pitch and roll as sine waves. yaw is 0
yaw = 0*t;
pitch = pi*sin(omega1*t);
roll =  pi*sin(omega2*t);

% For testing
%plot(t,pitch);
%hold on
%plot(t,roll);

%% Get Magnetometer and Star Tracker Readings Over 24h

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

%% Apply NLLS to Determine Attitude for each timestep
% Constants
tol = 0.01;
%Weights matrix same size as readings
weight_star = 1;
weight_mag = 0.01;

%find number of star readings, includes X,Y,Z
[num_times, c, num_star_readings] = size(Star_LGCV);

%*************************************** TEMP ****************************
num_times = 1;
%*************************************************************************

%loop for the number of times
for t = 1:num_times
    %init X vector and other variables
    X_vector = [0;0;0];    %yaw, pitch, roll
    delta_x = 10e9;
    itter = 0;
    weights = 0;
    % Perform NLLS
    while norm(delta_x) > tol
        %build jacobian for magnetometer readings
        k = 1;
        %calc estimated body frame coordinates using guess attitude
        y_est(3*k-2:3*k,1) = LGCV_to_Body(X_vector,Mag_LGCV(t,:));
        %build jacobian
        H(3*k-2:3*k,:) = Get_Jacobian_Attitude(X_vector, Mag_LGCV(t,:));
        %organise measurements
        y_meas(3*k-2:3*k,1) = Mag_Body(t,:);
        %build diagonals of weight matrix
        weights(3*k-2:3*k) = [weight_mag, weight_mag, weight_mag];
        
        %build jacobian for star tracker readings
        for k = 2:num_star_readings
            %calc estimated body frame coordinates using guess attitude
            y_est(3*k-2:3*k,1) = LGCV_to_Body(X_vector,Star_LGCV(t,:,k));
            %build jacobian
            H(3*k-2:3*k,:) = Get_Jacobian_Attitude(X_vector, Star_LGCV(t,:,k));
            %organise measurements
            y_meas(3*k-2:3*k,1) = Star_Body(t,:,k);
            %build diagonals of weight matrix
            weights(3*k-2:3*k) = [weight_star, weight_star, weight_star];
        end
        
        %get delta yo
        delta_y0 = y_meas - y_est;
        %make weights matrix
        W = diag(weights);
        %get delta x
        delta_x = (H'*W*H)\H'*W*delta_y0;
        %add to the X vector
        X_vector = X_vector + delta_x;
        
        %break condition for itterations
        itter = itter + 1;
        if itter > 10
            break;
        end
        
    end     %end tolerance loop for NLLS
    disp(X_vector)
end

%% Plots
figure
plot(t,pitch,'b')
hold on
plot(t,YPR_sun(:,2),'r');
hold on
plot(t,YPR_star(:,2),'k');
title('Pitch vs Time')
legend('Real Values','Sun Sensor','Star Tracker')