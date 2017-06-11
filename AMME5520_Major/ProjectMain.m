%%% Base file for AMME5520 Assignment 2
% Navigation and control through an obstacle field. The objective is to go
% from left to right.
clc
close all
clear all

%************************* PARAMETERS *************************************%
Width = 300;
Height = 150;

% Randomly generate some obstacles (number, average size as parameters)
numObst = 20; % Control the number of obstacles, e.g. 10, 20, 30, 40.
Adim = 15; % Control the "average" size of the obstacles.
buffer = 0.8;    %buffer zone in percent

% Starting/Ending Points in form [X,Y]
starting_point = [0,100];
ending_point = [250,100];

% Parameters for Path Planning
drawrealtime = 0;   % do PRM showing each step (looks cool!)
N = 200; %number of points/nodes required
K = 3;  %number of nearest nodes to connect to

% Parameters for Closed Loop Simulation
h = 0.01; % 100Hz sample time. Modify as desired.

%*************************************************************************%

%% Obstacle Generation
As = cell(numObst,1);
cs = cell(numObst,1);

dimensions = [0 Width 0 Height];

% Create obstacles only in the middle 60%
leftlim = 0.2*(dimensions(2)-dimensions(1));
rightlim = 0.8*(dimensions(2)-dimensions(1));

for k = 1:numObst
    % Generate ellipse in the form (x-c)'A(x-c)=1
    L = randn(2,2);
    As{k} = (0.4*eye(2)+ L'*L)/Adim^2;
    tmp = rand(2,1);
    cs{k} = [leftlim+(rightlim-leftlim)*tmp(1);dimensions(3)+(dimensions(4)-dimensions(3))*tmp(2)];

    %Add offset to ellipse, is proportional
    Ass(:,:,k) = As{k};
    Ass(:,:,k) = Ass(:,:,k)*buffer;        
end

%% Students add code here for planning

% - Students to modify ComputePath, to include collision checker, PRM, and
% path planner.

% [min_dist, path] = ComputePath(drawrealtime,N,K,Width,Height,dimensions,As,Ass,cs,starting_point,ending_point);
% save('workspace_11am')
 clear all
 load('workspace_11am')
 
 % Parameters for Kalman
C_kal = eye(8);

sigma_avg = 0.01;
sigma_X = 1/3;  % error from GPS
sigma_th = 0.1; % error from IMU/Gyro
sigma_thd = 0.01;   % error from IMU/Gyro
variance = [sigma_X; sigma_X; sigma_th; sigma_avg; sigma_avg; sigma_thd; sigma_avg; sigma_avg];

R = eye(8)*sigma_avg;
R(1,1) = sigma_X^2;
R(2,2) = sigma_X^2;
R(3,3) = sigma_th^2;
R(6,6) = sigma_thd^2;

P = eye(8)*1e-6;
%% Simulate Closed-loop system
% students to modify functions

h = 0.01; % 100Hz sample time. Modify as desired.
velocity = 5;   %m/s constant speed of the drone

x0 = zeros(8,1);  %Initial condition.
x0(1:2) = path(:,1);

stop = 0;
ts = 0;
xs = x0;

% This variable can be used to collect any parameters for simulation
m = 0.612;  %kg
g = 9.81;
L = 0.3;      %m length of copter
I = 3.03e-3;    %kgm^2
freq = 100; %Hz
delta_T = 1/freq;   % seconds

% Current State
xt = x0;

% Control Input
ut = [m*g/2; m*g/2];    %set initial values to hover at equillibrium

% Reference Control Input
u0 = [m*g/2; m*g/2];

% The parameters of the system
dynparams = [m,g,L,I,freq];

% Discretise Curve to Obtain desired state at each timestep.
%X_desired = equalspacing(path,delta_T  *velocity,velocity);
X_desired = Get_X_desired(path, delta_T*velocity, velocity, dynparams);
[~,num_steps] = size(X_desired);
k = 1;
while (stop ~= 1)
    % Assume K = Previous Step k == k - 1
    % K + 1 = Current Step we computing for k + 1 == k
    
    % Current state.
    xt = xs(:,k);
   
    ut = ComputeControl(u0, xt, X_desired(:,k), dynparams);
      
    % Use Runge Kutta to approximate the nonlinear dynamics over one time
    % step of length h.
    xs(:,k+1) = RungeKutta4(@QuadDynamics, xt, ut, 0, h, dynparams); 
    
    % Measurement of step K+1
    yt_k1 = meas(xs(:,k+1),variance);
    
    % Kalman filter state estimator 
    [P(:,:,k+1),xs(:,k)] = Kalman_Filter(xs(:,k), P(:,:,k), yt_k1, C_kal, R);
                   
    ts(k+1) = ts(k)+h;
    us(:,k) = ut;
    
    if (k>=num_steps)
        % if not reached, continue to try to reach?
        
        stop = 1;  % Stop after 100 time steps. You will need to change this.
        % Should stop after goal reached or collision.
    end
    
    k = k+1;
end

%% Removing outliers
for i = 1:num_steps-1
    if abs(xs(1,i)) > 300
        xs(:,i) = xs(:,i+1);
        disp(i)
    end
end

% Tracking Plot
figure
plot(xs(1,:),xs(2,:),'r')
hold on
plot(X_desired(1,:),X_desired(2,:),'b')
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
plot(path(1,:),path(2,:),'k.')
plot(path(1,1),path(2,1),'pk','MarkerSize',10)
plot(path(1,end),path(2,end),'pk','MarkerSize',10)
for k = 1:numObst    
    Ellipse_plot(As{k},cs{k},'b');
    Ellipse_plot(Ass(:,:,k),cs{k},'r');
    hold on   
end
title('Path Finding and Following')
xlabel('X axis (m)')
ylabel('Y axis (m)')

disp('done')