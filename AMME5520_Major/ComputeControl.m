function [u,params] = ComputeControl(y,params)
% Students to add code here. You can use the "params" variable to store
% internal controller states, if necessary.

% Get K controller using code from Assignment 1
% The model is linearised about the equillibirum
K = LQR_Get_K(params);

u = -K*y;

%u = [0;0];