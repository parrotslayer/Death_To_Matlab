function [ xdot ] = QuadDynamics(xt,ut,t, dynparams)
% Dynamics stub for AMME5520 Assignment 2. Students to fill in.
%clc 
%clear all
%%%%%%%%%%%% INPUTS ***********************
% m = 0.612;  %kg
% g = 9.81;
% L = 0.3;      %m length of copter
% I = 3.03e-3;    %kgm^2
% freq = 100; %Hz
% delta_T = 1/freq;   % seconds
% 
% % Current State
% xt = [10,10,10,0.5,0.5,0.5,10,10];
% 
% % Control Input
% ut = [m*g/2, m*g/2];
% 
% % Current Time
% t = 0;
% 
% % The parameters
% dynparams = [m,g,L,I,freq];
%***********************************************

% Give variables useful names
m = dynparams(1);
g = dynparams(2);
L = dynparams(3);
I = dynparams(4);
freq = dynparams(5);

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

T1_d_new = 2*r1 - 2*T1;
T2_d_new = 2*r2 - 2*T2;

% integrate T1_d to get T1 using 1st order Taylor Series
T1_new = T1 + delta_T*T1_d_new;
T2_new = T2 + delta_T*T2_d_new;

% Use new T1 and T2 to get Th''.
th_dd_new = 1/I*(L*(T1_new-T2_new));

% Get th' using new T1, T2, th''
th_d_new = th_d + delta_T*th_dd_new;

% Get th using new th'
th_new = th + delta_T*th_d_new;

% Use new T1 and T2 to get X'', Y''.
x_dd_new = 1/m*( (T1_new + T2_new)*sin(th_new) );
y_dd_new = 1/m*( -m*g + (T1_new + T2_new)*cos(th_new) );


% Integrate X'' to X' using simple 1st order Taylor Barker
x_d_new = x_d + delta_T*x_dd_new;
y_d_new = y_d + delta_T*y_dd_new;

% Output the derivative of the state vector
xdot = [x_d_new; y_d_new; th_d_new; x_dd_new; y_dd_new; th_dd_new; T1_d_new; T2_d_new];

end

