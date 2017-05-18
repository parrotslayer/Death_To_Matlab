% AERO4701 Space 3 Assignment 2 Question 1
clc
%close all
clear all

%% Part A
%Constants
rad = 6378137; %earth radius km
deg2rad = pi/180;
rad2deg = 180/pi;
G = 6.67408e-11;
M_earth = 5.972e24;
mu = G*M_earth;

% location of ground station
lat = -34.76; %S
long = 150.03; %E
height = 680;   %m

%% Import Data and initialise variables
GPS_ephem = dlmread('GPSsat_ephem.txt');
GPS_F2 = dlmread('GPS_pseudorange_F2.txt');
UAVPosition_F2 = dlmread('UAVPosition_F2.txt');

ephsat = GPS_ephem(:,2:8);
% Time of last vernal equinox
T_VE = 7347737.336;
% Find the last measurement of epoch, so we begin from there
lastepoch = max(ephsat(:,7));
% Check which observations from the pseudorange have same time
unique_times = unique(GPS_F2(:,1));
num_times = length(unique_times);

%% For a specific timestep
for i = 1 : num_times
    % Time for ECI calculations, since last satellite epoch measurement
    current_time = min(unique_times) + (i - 1);
    % Time since the last equinox, for ECEF coordinate calculation
    tsinceequinox = min(unique_times) + (i - 1) - T_VE;
    % Get the ECEF coordinates of all satellites at this specific time
    for n = 1 : 31
        timesinceepoch = current_time - ephsat(n,7);
        ECI_sat(:,n) = Orbit_to_ECI(ephsat(n,1),ephsat(n,2),...
            ephsat(n,3),ephsat(n,4),...
            ephsat(n,5),ephsat(n,6),...
            timesinceepoch);
        ECEF_sat(i,:,n) = ECI_to_ECEF([ECI_sat(:,n);tsinceequinox]);
    end
end

%% Organise GPS_F1 data
% Creates the GPS_pseudo matrix that stores the pseudorange readings in the first page
% and the corresponding source satellite in the 2nd page
[rows,cols] = size(GPS_F2);
sat_count = zeros(num_times,1);
% for each data set store GPS data
for i=1:rows
    %normalised time wrt t_0, starts at 1
    t = GPS_F2(i,1) - GPS_F2(1,1) + 1;
    %increment number of satellites seen
    sat_count(t) = sat_count(t) + 1;
    %store pseudo range
    GPS_pseudo(t,sat_count(t),1) = GPS_F2(i,3);
    %store satellite number on page 2
    GPS_pseudo(t,sat_count(t),2) = GPS_F2(i,2);
end

%% Nonlinear Least Squares
tolerance = 1;    %tolerance for convergence
itterations = 10;    %max number of itterations
    not_removed = 0;

[rows_GPS,cols_GPS,pages_GPS] = size(GPS_pseudo);
%do nonlinear regression for each timestep of the UAV
t=1;    %debuging
for t=1:num_times
    %***********************************************************************
    
    % Loop through all combinations where one satellite is removed
    for h = 1:sat_count(t)
                
        X_vector = [0;0;0;0];    %initial estimate
        %create 1s in the H matrix of according size
        clear H;
        H = ones((sat_count(t)-1),4);    %create (n-1)x4 matrix of 1s
        %create delta_p matrix
        clear delta_p
        delta_p = zeros(1);
        
        %calculate delta_x and repeat until converges or enough itterations
        for j=1:itterations
            f = 0;
            %% Find out which satellite is bad
            %calculate pseudorange and H for each satellite in view
            for k=1:sat_count(t)
                %make sure not to use measurement of satellite we are excluding
                if k ~= h
                    f = f+1;
                    %get corresponding satellite to measurement f
                    sat_k = GPS_pseudo(t,k,2);
                    % Get the ECEF coordinates for current satellite
                    x_sat = ECEF_sat(t,1,sat_k);
                    y_sat = ECEF_sat(t,2,sat_k);
                    z_sat = ECEF_sat(t,3,sat_k);
                    
                    % Get the difference in measured, and predicted UAV coordinates (X)
                    norm_x = x_sat - X_vector(1);
                    norm_y = y_sat - X_vector(2);
                    norm_z = z_sat - X_vector(3);
                    
                    % Calculate values for the jacobian, H
                    dpdx =  - norm_x /  sqrt( norm_x ^ 2 + norm_y ^ 2  + norm_z ^ 2);
                    dpdy =  - norm_y /  sqrt( norm_x ^ 2 + norm_y ^ 2  + norm_z ^ 2);
                    dpdz =  - norm_z /  sqrt( norm_x ^ 2 + norm_y ^ 2  + norm_z ^ 2);
                    
                    % Store values in row k of the Jacobian H
                    H(f,:) = [dpdx , dpdy, dpdz, 1];
                    
                    % get the estimated range for row k
                    range_est(f) = sqrt( (norm_x)^2 + (norm_y)^2 + (norm_z)^2 ) + X_vector(4);
                    
                    %calculate residuals for kth row = measured R - estimated R
                    delta_p(f,:) = GPS_pseudo(t,k,1)-range_est(f);
                end
            end
            
            %calculate the residuals after H and range_est are built
            delta_x = inv(transpose(H)*H) * transpose(H) * delta_p;
            
            %get the degree of precisions
        V = inv(transpose(H)*H);
        %check rank = 4?
        GDOP(t,h) = sqrt(V(1,1)+V(2,2)+V(3,3)+V(4,4));
        PDOP(t,h) = sqrt(V(1,1)+V(2,2)+V(3,3));
        HDOP(t,h) = sqrt(V(1,1)+V(2,2));
        VDOP(t,h) = sqrt(V(3,3));
        TDOP(t,h) = sqrt(V(4,4));
        
        %discard bad entries
        if rank(V) ~= 4
            GDOP(t,h) = NaN;
            PDOP(t,h) = NaN;
            HDOP(t,h) = NaN;
            VDOP(t,h) = NaN;
            TDOP(t,h) = NaN;
        end
            
        %check if within tolerance
            if delta_x < tolerance
                break;
            end
            
            %set new X as the residuals added on top
            X_vector = X_vector + delta_x;
        end   %end j / itterations loop
        %after getting the X_vector for a constellation without 'h'
        X_combinations(h,:) = X_vector;
        
    end     %end h loop
    
    %% We will look for ones with very bad clock biases
    thresh_CB = 1.5;   %threshold used to remove outliers for clock bias
    
    IQR_CB = iqr(X_combinations(:,4));
    Quantiles = quantile(X_combinations(:,4),[.25 .50 .75]);
    CB_upper = Quantiles(3) + thresh_CB*IQR_CB;
    CB_lower = Quantiles(1) - thresh_CB*IQR_CB;
    
    %discard values that are bad
    bad_count = 0;
    bad_sats = 0;
    testing = 0;
    for z=1:sat_count(t)
        if X_combinations(z,4) > CB_upper || X_combinations(z,4) < CB_lower
            %these are bad satellites so store them
            bad_count = bad_count + 1;
            %store the number of the bad satellite
            bad_sats(bad_count) = z;
            Worst_Sat = z;
        %if cant detect error, simply take the 1st one
        testing = 1;
        end     
    end %end for z 
    
    %if nothing bad is found
    if testing == 0
            Worst_Sat = 1;
            not_removed = not_removed + 1;
            plot_not_removed(not_removed) = t;
    end
    
    %if more than 2 found, find out which one is worse to remove
    if bad_count > 1
        min_error = 0;
        for p = 1:bad_count
            error = Quantiles(2) - X_combinations(bad_sats(p));
            %check if this is the largest error
            if error > min_error
                min_error = error;
                Worst_Sat = bad_sats(p);
            end
        end
    end
    
    % Now want to check which combination of satellite gives different
    % X_vector as this will be the one with the faulty satellite.
    X_vector_good = X_combinations(Worst_Sat,:);
    GDOP_plot = GDOP(t,Worst_Sat);
    PDOP_plot = PDOP(t,Worst_Sat);
    HDOP_plot = HDOP(t,Worst_Sat);
    VDOP_plot = VDOP(t,Worst_Sat);
    TDOP_plot = TDOP(t,Worst_Sat);

    %% store the final Delta_X and convert it to LLH,LGCV and RAE
    %convert UAV coordinates to LLH, LGCV, RAE
    UAV_ECEF_CB(:,t) = X_vector_good;
    
    %convert ECEF to LLH
    UAV_LLH(:,t) = ECEF_to_LLH(UAV_ECEF_CB(1,t), UAV_ECEF_CB(2,t), UAV_ECEF_CB(3,t));
    
    %convert ECEF to LGCV
    UAV_LGCV(:,t) = ECEF_to_LGCV(lat, long, height, UAV_ECEF_CB(1,t),...
        UAV_ECEF_CB(2,t), UAV_ECEF_CB(3,t));
    
    %convert LGCV to RAE
    UAV_RAE(:,t) = LGCV_to_RAE(UAV_LGCV(:,t));
end   %end t

% Remove bad data from the series using Clock Bias
thresh_CB = 20;   %threshold used to remove outliers for clock bias

IQR_CB = iqr(UAV_ECEF_CB(4,:));
Quantiles = quantile(UAV_ECEF_CB(4,:),[.25 .50 .75]);
CB_upper = Quantiles(3) + thresh_CB*IQR_CB;
CB_lower = Quantiles(1) - thresh_CB*IQR_CB;

%discard values that are bad
for t=1:num_times
    if UAV_ECEF_CB(4,t) > CB_upper || UAV_ECEF_CB(4,t) < CB_lower
        UAV_ECEF_CB(:,t) = NaN;
        UAV_LGCV(:,t) = NaN;
        UAV_RAE(:,t) = NaN;
        UAV_LLH(:,t) = NaN;
    end
end

%% Remove bad data from the series using LGCV Z values
thresh_Z = 10;   %threshold used to remove outliers for LGCV Z

IQR_Z = iqr(UAV_ECEF_CB(4,:));
Quantiles = quantile(UAV_LGCV(3,t),[.25 .50 .75]);
Z_upper = Quantiles(3) + thresh_Z*IQR_Z;
Z_lower = Quantiles(1) - thresh_Z*IQR_Z;

%discard values that are bad
for t=1:num_times
    if UAV_LGCV(3,t) > Z_upper || UAV_LGCV(3,t) < Z_lower
        UAV_ECEF_CB(:,t) = NaN;
        UAV_LGCV(:,t) = NaN;
        UAV_RAE(:,t) = NaN;
        UAV_LLH(:,t) = NaN;
    end
end

%% Plot the stuff
%close all
clc

%plot 3D overhead plot in LGCV
figure
plot3(UAV_LGCV(1,:),UAV_LGCV(2,:),UAV_LGCV(3,:),'.')
title('3D Overhead Plot in LGCV')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)') %errors in Z axis, lots of spread

hold on
plot3(UAVPosition_F2(:,2),UAVPosition_F2(:,3),UAVPosition_F2(:,4),'r')
legend('Estiamted','Real Position LGCV')

%plot 2D overhead plot
figure
for i = 1:num_times
    [plotx(i), ploty(i)] = overheadplotcoords(UAV_RAE(2,i), UAV_RAE(3,i));
end

plot(plotx,ploty,'.');
title('2D Overhead Plot ')
xlabel('Azimuth (degrees)')
ylabel('Elevation (degrees)')

%Plot LGCV Z vs time with reference data
figure
plot(abs(UAV_LGCV(3,:)));
title('Z vs Time')
xlabel('Time Step [t]')
ylabel('Elevation (m)')
hold on
plot(abs(UAVPosition_F2(:,4)));
hold on
plot(plot_not_removed,750,'k.')
legend('Estiamted Elevation','Real Position LGCV','No Bad Satellite Found')


% Plot Clock bias and number of satellites
figure
ax1 = subplot(3,1,1);
plot(ax1,unique_times, UAV_ECEF_CB(4,:));
title(ax1,'Clock Bias vs Time')
xlabel(ax1,'Epoch Time (seconds)')
ylabel(ax1,'Clock Bias (m)')

ax2 = subplot(3,1,2);
plot(ax2,unique_times, sat_count);
title(ax2,'Number satellites visible vs Time')
xlabel(ax2,'Epoch Time (seconds)')
ylabel(ax2,'Number of Visible Satellites')

%plot altitude vs time
ax3 = subplot(3,1,3);
plot(unique_times, UAV_LLH(3,:));
title(ax3,'Elevation vs Time')
xlabel(ax3,'Epoch Time (seconds)')
ylabel(ax3,'Elevation (m)')

%plot DOPs
figure
w = Worst_Sat;
plot(unique_times,GDOP(:,w));
hold on
plot(unique_times,PDOP(:,w));
hold on
plot(unique_times,HDOP(:,w));
hold on
plot(unique_times,VDOP(:,w));
hold on
plot(unique_times,TDOP(:,w));
hold on
title('Dilution of Precision vs Time')
xlabel('Epoch Time (seconds)')
ylabel('Dilution of Precision')
legend('GDOP','PDOP','HDOP','VDOP','TDOP')
%% Get RAE of 1st and last timestep for all visible satellites
% get RAE for 1st timestep
t=1;
for k=1:sat_count(t)
    %get corresponding satellite to measurement k
    sat_k = GPS_pseudo(t,k,2);
    % Get the ECEF coordinates for current satellite
    x_sat = ECEF_sat(t,1,sat_k);
    y_sat = ECEF_sat(t,2,sat_k);
    z_sat = ECEF_sat(t,3,sat_k);
    
    ECEF_1 = ECEF_sat(t,:,sat_k);
    LGCV_1 = ECEF_to_LGCV(lat,long,height,ECEF_1(1),ECEF_1(2),ECEF_1(3));
    RAE_1(:,k) = LGCV_to_RAE(LGCV_1);
end

% get RAE for last timestep
t=num_times;
for k=1:sat_count(t)
    %get corresponding satellite to measurement k
    sat_k = GPS_pseudo(t,k,2);
    % Get the ECEF coordinates for current satellite
    x_sat = ECEF_sat(t,1,sat_k);
    y_sat = ECEF_sat(t,2,sat_k);
    z_sat = ECEF_sat(t,3,sat_k);
    
    ECEF_2 = ECEF_sat(t,:,sat_k);
    LGCV_2 = ECEF_to_LGCV(lat,long,height,ECEF_2(1),ECEF_2(2),ECEF_2(3));
    RAE_2(:,k) = LGCV_to_RAE(LGCV_2);
end

%% Plot Azimuth and Elevation of the Satellites visible
% plotting of an azimuth-elevation overhead plot. Also use
% "overheadplotcoords" to transform azimuth and elevation angles into plot
% coordinates for this figure

% Setup Figure
figure
hold on

% Plot Markings on Figure
tincs = 0:0.2:(2*pi + 0.2);
xpol(1,:) = 30*sin(tincs);
xpol(2,:) = 60*sin(tincs);
xpol(3,:) = 90*sin(tincs);
ypol(1,:) = 30*cos(tincs);
ypol(2,:) = 60*cos(tincs);
ypol(3,:) = 90*cos(tincs);
plot([0,0],[-90,90])
plot([-90,90],[0,0])
plot(xpol(1,:),ypol(1,:))
plot(xpol(2,:),ypol(2,:))
plot(xpol(3,:),ypol(3,:))

% Plot Directions and Numbers on Figures
text(0,90,'North')
text(-90,0,'East')
text(0,-90,'South')
text(90,0,'West')
text(0,0,'90^o')
text(30*0.7,30*0.7,'60^o')
text(60*0.7,60*0.7,'30^o')
text(90*0.7,90*0.7,'0^o')


plotx = 0;
ploty = 0;

for i=1:sat_count(1)
    [plotx(i), ploty(i)] = overheadplotcoords(RAE_1(2,i), RAE_1(3,i));
end
plot(plotx,ploty,'o');

hold on
for i=1:sat_count(num_times)
    [plotx(i), ploty(i)] = overheadplotcoords(RAE_2(2,i), RAE_2(3,i));
end
plot(plotx,ploty,'o');
title('Azimuth Elevation Plot of Visible Satellites at Start (Blue) and End (Red)')
xlabel('Azimuth (Degrees)')
ylabel('Elevation (Degrees)')