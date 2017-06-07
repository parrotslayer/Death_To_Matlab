%function [ xdot ] = QuadDynamics(xt,ut,t, dynparams)
% Dynamics stub for AMME5520 Assignment 2. Students to fill in.
clc 
clear all
%%%%%%%%%%%% INPUTS ***********************
m = 0.612;  %kg
g = 9.81;
L = 0.3;      %m length of copter
I = 3.03e-3;    %kgm^2
freq = 100; %Hz
delta_T = 1/freq;   % seconds

% Current State
xt = [10,10,10,0.5,0.5,0.5,10,10];

% Control Input
ut = [m*g/2, m*g/2];

% Current Time
t = 0;

% The parameters
dynparams = [m,g,L,I,freq];
%***********************************************

% Give variables useful names
m = dynparams(1);
g = dynparams(2);
L = dynparams(3);
I = dynparams(4);
r1 = ut(1);
r2 = ut(2);
x = xt(1);
y = xt(2);
th = xt(3);
x_d = xt(4);
y_d = xt(5);
th_d = xt(6);
T1 = xt(7);
T2 = xt(8);

delta_T = 1/freq;   % seconds

% Get K controller using code from Assignment 1
% The model is linearised about the equillibirum
K = LQR_Get_K();

disp(K)

x_dd = 1/m*( (T1+T2)*sin(th) );
y_dd = 1/m*( -m*g + (T1+T2)*cos(th) );
th_dd = 1/I*(L*(T1-T2));
T1_d = 2*r1 - 2*T1;
T2_d = 2*r2 - 2*T2;

% Integrate X'' to X' using simple 1st order
x_d = x_d + delta_T*x_dd;
y_d = y_d + delta_T*y_dd;

% Output the derivative of the state vector
xdot = [x_d, y_d, th_d, x_dd, y_dd, th_dd, T1_d, T2_d];

%end

