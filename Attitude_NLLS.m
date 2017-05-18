% NLLS of attitude
clc
close all
clear all

timeperiod = 24*60*60; %1day
global step
step = 100;         %time step for looping calculations

%% Get Estimated and Real Satellite positions over 24h
orbit_params(1)= 7162*1000;               %a - semi major axis meters
orbit_params(2)=0.0000872;                %e - eccentricity deg
orbit_params(3)=98.7401;                  %inc - inclination degrees
orbit_params(4)=142.1145;                 %Omega - degrees
orbit_params(5)=25.5213;                  %omega - degrees
orbit_params(6)=283.570800000000;         %Mo - Mean Anomaly
orbit_params(7)=2452716.5000;             %Julian day = 2452716.5000 [days] 

% Use code from Assignment 2 to model estimated orbit and real orbit
[ECI_sat,ECI_est,LGCV_est] = Get_Estimated_ECEF(orbit_params);

%% Simulate the Attitude with Noise
t = 1:timeperiod;   %generate times
omega1 = 0.001;     %freq of pitch
omega2 = 0.005;     %freq of yaw

%model pitch and roll as sine waves. yaw is 0
yaw = 0*t;
pitch = pi*sin(omega1*t);
roll =  pi*sin(omega2*t);

% For testing
plot(t,pitch);
hold on
plot(t,roll);

%% Get Magnetometer and Star Tracker Readings Over 24h
for t = 1:step:time_period
    
    
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