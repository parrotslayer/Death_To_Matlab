% clc
% clear all

function [P,x_hat] = Kalman_Filter(x_hat_minus, P_prev, yt, C, R)

% **************************** Constants **********************************
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
Q = eye(8);
% **********************************************************************

% % Prediction Stage Estimates x^-(k)
% x_hat_minus = A*x_hat_prev + B*ut;
P_minus = A*P_prev*A' + Q

% Equation 1. Calculate Kalman Gain K
K = (P_minus*C')\( C*P_minus*C' + R );
% K(isnan(K)) = 0;

% Equation 2. Calculate estimated x_hat
x_hat = x_hat_minus + K*(yt - C*x_hat_minus);

% Equation 3. Calculate new Covariance
P = (eye(8) - K*C)*P_minus;

disp('K = ')
disp(K)
%  disp('X hat')
%  disp(x_hat)
%  disp('P')
%  disp(P)
%  disp('done')
end
