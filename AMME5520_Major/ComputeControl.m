% Takes in the previous control (Y), current state and desired current
% state and paramters.
% Outputs the required u (y) to control the drone.
function [u,params] = ComputeControl(u0, X_current, X_desired, params)

% Get K controller using code from Assignment 1
% The model is linearised about the equillibirum
K = LQR_Get_K(params);

delta_X = X_current - X_desired;

u = u0 - K*delta_X;

