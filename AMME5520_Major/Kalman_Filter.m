clc
clear all

%function [P,x_hat] = Kalman_Filter(x_hat_prev, P_prev, yt, ut, C_kal, R)
%*********************** Constants *******************************
m = 0.612;  %kg
g = 9.81;
L = 0.3;      %m length of copter
I = 3.03e-3;    %kgm^2

A = [0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 g 0 0 0 0 0;
    0 0 0 0 0 0 1/m 1/m;
    0 0 0 0 0 0 L/I -L/I;
    0 0 0 0 0 0 -2 0;
    0 0 0 0 0 0 0 -2];

% Defining the B matrix
B = [0 0;
    0 0;
    0 0;
    0 0;
    0 0;
    0 0;
    2 0;
    0 2];

% **************************** INPUTS **********************************
C_kal = zeros(8);
C_kal(1,1) = 1; % X
C_kal(2,2) = 1; % Y
C_kal(3,3) = 1; % th
C_kal(6,6) = 1; % th;

R = zeros(8);
sigma_X = 0.01;
sigma_Y = 0.01;
sigma_Th = 0.01;
sigma_Th_dot = 0.01;

R(1,1) = sigma_X;
R(2,2) = sigma_Y;
R(3,3) = sigma_Th;
R(6,6) = sigma_Th_dot;

x_hat_prev = [10;10;0;0;0;0;3;3];
P_prev = eye(8)*1e-1;
y = 1; 
ut = [3;3];
% **********************************************************************

% Prediction Stage Estimates x^-(k)
x_hat_minus = A*x_hat_prev + B*ut;

% Equation 1. Calculate Kalman Gain K
K = (P_prev*C_kal')/( C_kal*P_prev*C_kal' + R );

% Equation 2. Calculate estimated x_hat
x_hat = x_hat_minus + K*(yt - C_kal*x_hat_minus);

% Equation 3. Calculate new Covariance
P = (eye(8) - K*C_kal)*P_prev;

disp(x_hat)
disp(P)
%end