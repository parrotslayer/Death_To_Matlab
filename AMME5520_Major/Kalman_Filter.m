% clc
% clear all

function [P,x_hat] = Kalman_Filter(x_hat_minus, P_prev, yt, ut, C, R)

% **************************** INPUTS **********************************
% C = zeros(8);
% C(1,1) = 1; % X
% C(2,2) = 1; % Y
% C(3,3) = 1; % th
% C(6,6) = 1; % th;
% C = eye(8);
% 
% R = zeros(8);
% sigma_X = 0.2;
% sigma_Y = 0.2;
% sigma_Th = 0.2;
% sigma_Th_dot = 0.2;
% noise = sigma_X;
% 
% R(1,1) = sigma_X^2;
% R(2,2) = sigma_Y^2;
% R(3,3) = sigma_Th^2;
% R(6,6) = sigma_Th_dot^2;
% 
% 
% P_prev = eye(8)*1e-10;
% 
% % mg/2
% ut = [3;3];
% ut = [12.6195594178204;1.37064440905142];
% 
% % From Project Main
% xs =[0;100;0;0;0;0;0;0];
% x_hat_minus = xs;
% yt = [0;100;0;0;0;0;0;0];
% yt = normrnd(yt,noise); 
% **********************************************************************

% % Prediction Stage Estimates x^-(k)
% x_hat_minus = A*x_hat_prev + B*ut;

% Equation 1. Calculate Kalman Gain K
K = (P_prev*C')\( C*P_prev*C' + R );
% K(isnan(K)) = 0;

% Equation 2. Calculate estimated x_hat
x_hat = x_hat_minus + K*(yt - C*x_hat_minus);

% Equation 3. Calculate new Covariance
P = (eye(8) - K*C)*P_prev;

 disp('X hat')
 disp(x_hat)
 disp('P')
 disp(P)
 disp('done')
%end