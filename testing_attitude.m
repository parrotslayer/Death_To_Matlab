%testing conversion matrix
clc
%close all
clear all

%% Testing Conversion Matrix
Star_Body(1,:,1) = [-0.0879; 0.5242; -0.6383];    %meas1
Star_Body(1,:,2) = [-0.3319; 0.3281; 0.2055];     %meas2
Star_Body(1,:,3) = [0.8465; 0.0540; 0.6485];      %meas3

Star_LGCV(1,:,1) = [0.2; 0.7; -0.4];   %meas 1
Star_LGCV(1,:,2) = [0.1; 0.3; 0.4];   %meas 1
Star_LGCV(1,:,3) = [0.7; -0.8; 0.1];   %meas 1

%magnetometers only take 1 reading
Mag_Body(1,:) = [-0.0879; 0.5242; -0.6383];      %meas1
%model magnetometer with errors
Error1 = [0.01, -0.02, 0.04];
Mag_Body(1,:) = Mag_Body(:,1) + Error1;
Mag_LGCV(1,:) = [0.2; 0.7; -0.4];   %meas 1

%t = 2
Star_Body(2,:,1) = [-0.0879; 0.5242; -0.6383];    %meas1
Star_Body(2,:,2) = [-0.3319; 0.3281; 0.2055];     %meas2
Star_Body(2,:,3) = [0.8465; 0.0540; 0.6485];      %meas3

Star_LGCV(2,:,1) = [0.2; 0.7; -0.4];   %meas 1
Star_LGCV(2,:,2) = [0.1; 0.3; 0.4];   %meas 1
Star_LGCV(2,:,3) = [0.7; -0.8; 0.1];   %meas 1

%magnetometers only take 1 reading
Mag_Body(2,:) = [-0.0879; 0.5242; -0.6383];      %meas1
%model magnetometer with errors
Error1 = [0.01, -0.02, 0.04];
Mag_Body(2,:) = Mag_Body(2,:) + Error1;
Mag_LGCV(2,:) = [0.2; 0.7; -0.4];   %meas 1


% t = 3;
%magnetometers only take 1 reading
Mag_Body(3,:) = [-0.0879; 0.5242; -0.6383];      %meas1
%model magnetometer with errors
Error1 = [0.01, -0.02, 0.04];
Mag_Body(3,:) = Mag_Body(3,:) + Error1;
Mag_LGCV(3,:) = [0.2; 0.7; -0.4];   %meas 1
Star_Body(3,:,:) = 0;
Star_LGCV(3,:,:) = 0;

% t = 4
Star_Body(4,:,1) = [-0.0879; 0.5242; -0.6383];    %meas1
Star_Body(4,:,2) = [-0.3319; 0.3281; 0.2055];     %meas2

Star_LGCV(4,:,1) = [0.2; 0.7; -0.4];   %meas 1
Star_LGCV(4,:,2) = [0.1; 0.3; 0.4];   %meas 1

%magnetometers only take 1 reading
Mag_Body(4,:) = [-0.0879; 0.5242; -0.6383];      %meas1
%model magnetometer with errors
Error1 = [0.01, -0.02, 0.04];
Mag_Body(4,:) = Mag_Body(4,:) + Error1;
Mag_LGCV(4,:) = [0.2; 0.7; -0.4];   %meas 1

deg2rad = pi/180;
True_Euler = [11*deg2rad; 32*deg2rad; -45*deg2rad];

%for testing
output = LGCV_to_Body(True_Euler,Star_LGCV(1,:,1));

%% Apply NLLS
% Constants
tol = 0.01;
%Weights matrix same size as readings
weight_star = 1;
weight_mag = 0.01;
%find number of star readings, includes X,Y,Z
[num_times, c, num_star_readings] = size(Star_LGCV);

%loop for the number of times
for t = 1:num_times
    %init X vector and other variables
    X_vector = [0;0;0];    %yaw, pitch, roll
    delta_x = 10e9;
    itter = 0;
    weights = 0;
        
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