%testing conversion matrix
clc
close all
clear all

%% Testing Conversion Matrix
Star_Body(1,:,1) = [-0.0879; 0.5242; -0.6383];    %meas1
Star_Body(1,:,2) = [-0.3319; 0.3281; 0.2055];     %meas2
Star_Body(1,:,3) = [0.8465; 0.0540; 0.6485];      %meas3

%add errors
Error1 = [0.01, -0.02, 0.04];
Star_Body(1,:,1) = Star_Body(1,:,1) + Error1;

Star_LGCV(1,:,1) = [0.2; 0.7; -0.4];   %meas 1
Star_LGCV(1,:,2) = [0.1; 0.3; 0.4];   %meas 1
Star_LGCV(1,:,3) = [0.7; -0.8; 0.1];   %meas 1

deg2rad = pi/180;
True_Euler = [11*deg2rad; 32*deg2rad; -45*deg2rad];

%for testing
output = LGCV_to_Body(True_Euler,Star_LGCV(1,:,1));

lgcv = Body_to_LGCV(True_Euler,Star_Body(1,:,2))
%% Apply NLLS
X_vector = [0;0;0];    %yaw, pitch, roll
tol = 0.0001;
itter = 0;
delta_x = 10e9;

t = 1;

%find number of readings, includes X,Y,Z
[r, c, num_star_readings] = size(Star_LGCV);

%Weights matrix same size as readings
weights = [1 1 1 1 1 1 1 1 1];
W = diag(weights);

while norm(delta_x) > tol
    for k = 1:num_star_readings
        %calc estimated body frame coordinates using guess attitude
        y_est(3*k-2:3*k,1) = LGCV_to_Body(X_vector,Star_LGCV(t,:,k));
        %build jacobian
        H(3*k-2:3*k,:) = Get_Jacobian_Attitude(X_vector, Star_LGCV(t,:,k));
        %organise measurements
        y_meas(3*k-2:3*k,1) = Star_Body(t,:,k);
    end
    
    %get delta yo
    delta_y0 = y_meas - y_est;
    
    %get delta x
    delta_x = (H'*W*H)\H'*W*delta_y0;

    X_vector = X_vector + delta_x
    
    itter = itter + 1;
    if itter > 10
       break;
    end
   
end
kek = X_vector*180/pi