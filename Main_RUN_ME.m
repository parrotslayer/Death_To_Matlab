% NLLS of attitude
clc
close all
clear all

time_period = 24*60*60; %1day
global step
global gs_num;  %the chosen ground station  
deg2rad = pi/180;
rad2deg = 180/pi;
radius_earth = 6378137; %from lecture slides
% Northern: Monday, 20 March 2017 at 10:29 UTC
julian_date17 = juliandate(datetime([2017,3,20,10,29,0]))*86400;

gs_names = {'Adelaide, South Australia','Hermitage, UK','Cape Caneveral, Florida'};
gs_llh = [-33.9284989  ,  138.6007456   , 45;...
    51.4520259, -1.2754 , 110  ;...
    28.3922182,  -80.6077132 , 2];

%% Get Estimated and Real Satellite positions over 24h
orbit_params(1)= 7162*1000;              %a - semi major axis meters    
orbit_params(2)=0.0000872;                %e - eccentricity deg
orbit_params(3)=98.7401;                  %inc - inclination degrees
orbit_params(4)=142.1145;                 %Omega - degrees
orbit_params(5)=25.5213;                  %omega - degrees
orbit_params(6)=283.570800000000;         %Mo - Mean Anomaly
orbit_params(7)=2457866.50000000;         %Julian Day (Epoch) Sunday 23/4/17 UT1 00:00:00

% Use code from Assignment 2 to model estimated orbit and real orbit
[Sat_ECI_true,Sat_ECEF_true,Sat_ECI_est, Sat_ECEF_est] = Orbit_Determination(orbit_params);
save('Orbit Determination_Data')

%load('Orbit Determination_Data')


%% Simulate the Attitude
step = 100;
t = 1:time_period;   %generate times
omega1 = 0.0005;     %freq of pitch
omega2 = 0.001;     %freq of yaw

%model pitch and roll as sine waves. yaw is 0
amplitude_attitude = 10;    %degrees
yaw = 0*t;
pitch = amplitude_attitude*deg2rad*sin(omega1*t);
roll =  amplitude_attitude*deg2rad*sin(omega2*t);

Attitude_Real = [roll;pitch;yaw];

%% Generate constellation of stars

%Greenwich sidereal time at UT[hours] 205.3166 [deg]
th_g0 = 205.3166*deg2rad; %greenwhich sidereal time at epoch
time_epoch = orbit_params(7) * 86400;  %time of epoch in seconds

%constants for generating constellation
radius = 10e10;      %some big number 
inc = orbit_params(3);
Omega = orbit_params(4);
omega = orbit_params(5);
num_stars = 300;
FOV = 60*deg2rad;

% Generate the constellation of stars
Star_Constellation_ECI = Generate_Random_Stars(num_stars,radius);

%scale down distance of the stars for plotting
scale = 1e3;

% Plot the ECI orbit and the star constellation
figure
plot3(Sat_ECI_est(1,:),Sat_ECI_est(2,:),Sat_ECI_est(3,:),'r.')
hold on
plot3(Star_Constellation_ECI(:,1)/scale,Star_Constellation_ECI(:,2)/scale,Star_Constellation_ECI(:,3)/scale,'.')
axis equal
title('ECI Plot of Orbit (RED) and Scaled Down Stars (BLUE)')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

%% Get LGCV and Body frame coordinates
% Errors
sigma_star = 0.001;      % direction
sigma_mag = 0.01;        % direction

%Preallocate for speed
Sat_LLH_est = zeros(time_period,3);
Mag_ECI = zeros(3,time_period);
Mag_ECEF = zeros(time_period,3);
Mag_LGCV = zeros(time_period,3);
Mag_RAE = zeros(time_period,3);
Mag_Body = zeros(time_period,3);

%24 is an arbitrary number of max stars seen
Star_ECI = zeros(24,3,time_period);
Star_ECEF = zeros(24,3,time_period);
Star_LGCV = zeros(24,3,time_period);
Star_RAE = zeros(24,3,time_period);
Star_Body = zeros(24,3,time_period);

for t = 1:step:time_period
%% Get Magnetometer Body and LGCV Readings with Noise
% Convert satellite's estimated ECEF position to est LLH
Sat_LLH_est(t,:) = ECEF_to_LLH(Sat_ECEF_est(1,t),Sat_ECEF_est(2,t),Sat_ECEF_est(3,t));

% Get unit vector of magnetometer readings
Mag_ECI(:,t) = Get_Mag(t,th_g0, Sat_ECI_est(:,t));

% Convert ECI to ECEF
current_time = time_epoch + t;    %time since last epoch for time = n
t_since_equinox = current_time - julian_date17;
Mag_ECEF(t,:) = ECI_to_ECEF([Mag_ECI(:,t); t_since_equinox]);

% Convert ECEF to LGCV wrt Satellite's position
Mag_LGCV(t,:) =  ECEF_to_LGCV(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
    Mag_ECEF(t,1), Mag_ECEF(t,2), Mag_ECEF(t,3));  

% Convert LGCV to RAE
Mag_RAE(t,:) = LGCV_to_RAE(Mag_LGCV(t,:));

% Add noise to RAE (Use rand because is faster than normrnd Gaussian Noise)
Mag_RAE(t,1) = Mag_RAE(t,1); 
Mag_RAE(t,2) = normrnd(Mag_RAE(t,2), sigma_mag); 
Mag_RAE(t,3) = normrnd(Mag_RAE(t,3), sigma_mag); 

% Convert RAE with noise back to LGCV
Mag_LGCV(t,:) = RAE_to_LGCV(Mag_RAE);

% Convert LGCV with noise to Body using Real angles
Mag_Body(t,:) = LGCV_to_Body(transpose(Attitude_Real(:,t)),Mag_LGCV(t,:));

%% Get Star Tracker Body and LGCV Readings with Noise
% Get unit vectors from the star trackers
temp = Get_Star_Tracker(Sat_LLH_est(t,:),t_since_equinox,Star_Constellation_ECI,FOV);
% Size of matrix changes so this code makes sure the matrix can be added
[r,c] = size(temp);
Star_LGCV(1:r,:,t) = temp; 

for k = 1:r
% Convert LGCV to RAE
Star_RAE(k,:,t) = LGCV_to_RAE(Star_LGCV(k,:,t));

% Add noise to RAE (Use rand because is faster than normrnd Gaussian Noise)
Star_RAE(k,1,t) = Star_RAE(k,1,t); 
Star_RAE(k,2,t) = normrnd(Mag_RAE(t,2), sigma_mag); 
Star_RAE(k,3,t) = normrnd(Mag_RAE(t,3), sigma_mag); 

% Convert RAE with noise back to LGCV
Star_LGCV(k,:,t) = RAE_to_LGCV(Star_RAE(k,:,t));

% Convert LGCV to Body using Real Attitudes
Star_Body(k,:,t) = LGCV_to_Body(Attitude_Real(:,t),Star_LGCV(k,:,t));     
end

% Progress Report
if t == 10001 || t==20001 ||t==30001||t==40001||t==50001||t==60001||t==70001||t==80001
    disp('Frame Conversion Timestep = ')
    disp(t)
end

end

%% Apply NLLS to Determine Attitude for each timestep
% Constants
tol = 0.01;
%Weights matrix same size as readings
weight_star = 1;
weight_mag = 1;
%find number of star readings, includes X,Y,Z
[r,c,num_times] = size(Star_LGCV);

%Preallocate for speed
Attitude_Est = zeros(num_times,3);
num_star_readings = zeros(num_times,1);
DOP_Roll = zeros(num_times,1);
DOP_Pitch = zeros(num_times,1);
DOP_Yaw = zeros(num_times,1);

%*************************************** TEMP ****************************
%num_times = 10001;
%*************************************************************************

%loop for the number of times
for t = 1:step:num_times
%    for t = 7201:step:7301

    %init X vector and other variables
    X_vector = [0;0;0];    %yaw, pitch, roll
    delta_x = 10e9;
    itter = 0;

    %clear matricies in case jacobian is a different size
    clear H
    clear y_est
    clear y_meas
    clear weights
    clear delta_y0
    
    % Get number of star readings for this timestep
    num_star_readings(t) = length(find(Star_LGCV(:,1,t)));

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
        for k = 2:num_star_readings(t)
            %calc estimated body frame coordinates using guess attitude
            y_est(3*k-2:3*k,1) = LGCV_to_Body(X_vector,Star_LGCV(k,:,t));
            %build jacobian
            H(3*k-2:3*k,:) = Get_Jacobian_Attitude(X_vector, Star_LGCV(k,:,t));
            %organise measurements
            y_meas(3*k-2:3*k,1) = Star_Body(k,:,t);
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
        
        % Check Rank of H
        if rank(H) ~= 3
            disp('Error, H is not Full Rank at Timestep = ');
            disp(t);
        end
        
        %Calculate DOP
        V = inv(H'*H);
        DOP_Roll(t) = sqrt(V(1,1));
        DOP_Pitch(t) = sqrt(V(2,2));
        DOP_Yaw(t) = sqrt(V(3,3));
        
        %break condition for itterations
        itter = itter + 1;
        if itter > 10
            break;
        end
        
    end     %end tolerance loop for NLLS
    Attitude_Est(t,:) = X_vector;
    
    if t == 10001 || t==20001 ||t==30001||t==40001||t==50001||t==60001||t==70001||t==80001
        disp('Attitude Determination Timestep = ')
        disp(t)
    end
end

%% Calculate Nadir Point and North/East Swath Edge
% Tracking Errors
sigma_tracking_errors = 100;    %meters

Sat_LLH_true = NaN(time_period,3);
Nadir_ECEF_true = NaN(time_period,3);
Nadir_ECEF_est = NaN(time_period,3);
Nadir_Error = NaN(time_period,3);
Nadir_Error_Magnitude = NaN(time_period,1);

for t = 1:step:time_period
%% Get True Position of the Nadir Point and Swath Edge
% Convert satellite's estimated ECEF position to est LLH
Sat_LLH_true(t,:) = ECEF_to_LLH(Sat_ECEF_true(1,t),Sat_ECEF_true(2,t),Sat_ECEF_true(3,t));

% Set magnitude of ECI vector to Radius of the Earth
Feature_ECI_true = Sat_ECI_true(:,t)*radius_earth/norm(Sat_ECI_true(:,t));

% Convert ECI to ECEF
current_time = time_epoch + t;    %time since last epoch for time = n
t_since_equinox = current_time - julian_date17;
Nadir_ECEF_true(t,:) = ECI_to_ECEF([Feature_ECI_true; t_since_equinox]);

% Convert ECEF to LGCV

% Convert LGCV to Body

% Rotate Body Frame by FOV/2 in the North Direction
th = FOV/2;
C_rot = [cos(th) 0 -sin(th);
    0 1 0;
    sin(th) 0 cos(th)];
% Body to LGCV

% LGCV to ECEF

%% Get Estimated Nadir Point
% Set magnitude of ECI vector to Radius of the Earth
Feature_ECI_est = Sat_ECI_est(:,t)*radius_earth/norm(Sat_ECI_est(:,t));

% Convert ECI to ECEF
Feature_ECEF_est = ECI_to_ECEF([Feature_ECI_est; t_since_equinox]);

% Convert ECEF to LGCV wrt Satellite's position
Feature_LGCV_est =  ECEF_to_LGCV(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
    Feature_ECEF_est(1), Feature_ECEF_est(2), Feature_ECEF_est(3));  

% Convert LGCV to Body using estimated angles
Feature_Body_est = LGCV_to_Body(transpose(Attitude_Est(t,:)),Feature_LGCV_est');

% Add Errors to Body readings
Nadir_Body_est = normrnd(Feature_Body_est,sigma_tracking_errors);

% Convert Body to LGCV
Nadir_LGCV_est = Body_to_LGCV(Attitude_Est(t,:),Nadir_Body_est);

% Convert LGCV to ECEF
Nadir_ECEF_est(t,:) = LGCV_to_ECEF(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
    Nadir_LGCV_est(1), Nadir_LGCV_est(2), Nadir_LGCV_est(3));

% Calculate Nadir Errors
Nadir_Error(t,:) = Nadir_ECEF_true(t,:) - Nadir_ECEF_est(t,:);
Nadir_Error_Magnitude(t) = norm(Nadir_Error(t,:));
end

%% Tracking Plots

figure
plot(1:time_period,Nadir_ECEF_true(:,1),'b');
hold on
plot(1:time_period,Nadir_ECEF_est(:,1),'r.');
title('X Nadir Est vs Real')

figure
plot(1:time_period,Nadir_Error_Magnitude,'k.')
title('Error Magnitude vs Time')

figure
subplot(2,1,1)
plot(Nadir_Error(:,1),Nadir_Error(:,2),'.');
title('Nadir Mapping Errors')
subplot(2,1,2)
plot(1:time_period,Nadir_Error(:,3),'.')
title('Down Error')

%% Generate Plots of Attitude
% Get rid of zero terms (no data because we skip steps)
%temp = num_star_readings(time_period);
num_star_readings(num_star_readings == 0) = NaN;
%num_star_readings(time_period) = temp;
Attitude_Est(Attitude_Est == 0) = NaN;

figure
ax1 = subplot(2,1,1);
plot(ax1,1:time_period,Attitude_Real(1,:),'b')
hold on
plot(ax1,1:time_period,Attitude_Est(1:time_period,1),'r.')
title(ax1,'Unfilted Roll')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Angle (Radians)')
ax2 = subplot(2,1,2);
plot(ax2,abs(num_star_readings),'k.');
title(ax2,'Number Visible Stars')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Number of Stars')

figure
subplot(2,1,1)
plot(1:time_period,Attitude_Real(2,:),'b')
hold on
plot(1:time_period,Attitude_Est(1:time_period,2),'r.')
title('Unfilted Pitch')
xlabel('Time (seconds)')
ylabel('Angle (Radians)')
subplot(2,1,2)
plot(abs(num_star_readings),'k.');
title('Number Visible Stars')
xlabel('Time (seconds)')
ylabel('Number of Stars')

figure
subplot(2,1,1)
plot(1:time_period,Attitude_Real(3,:),'b')
hold on
plot(1:time_period,Attitude_Est(1:time_period,3),'r.')
title('Unfilted Yaw')
xlabel('Time (seconds)')
ylabel('Angle (Radians)')
subplot(2,1,2)
plot(abs(num_star_readings),'k.');
title('Number Visible Stars')
xlabel('Time (seconds)')
ylabel('Number of Satellites')

figure
subplot(4,1,1)
plot(1:time_period,DOP_Roll)
title('DOP Roll')
subplot(4,1,2)
plot(1:time_period,DOP_Pitch)
title('DOP Pitch')
subplot(4,1,3)
plot(1:time_period,DOP_Yaw)
title('DOP Yaw')
subplot(4,1,4)
plot(1:time_period,num_star_readings,'k.');
title('Number of Visible Stars')
xlabel('Time (Seconds)')