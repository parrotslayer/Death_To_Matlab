function [ output ] = QuadDynamics(xt,ut,t, dynparams)
% Dynamics stub for AMME5520 Assignment 2. Students to fill in.
m = 0.612;
g = 9.81;
I = 3.03e-3;
L = 0.3;
% x        = xt(1);
% y        = xt(2);
theta    = xt(3);
xdot     = xt(4);
ydot     = xt(5);
thetadot = xt(6);
T1       = xt(7);
T2       = xt(8);
% Thrust Calcs
T1dot = 2 * ut(1) - 2 * T1;
T2dot = 2 * ut(2) - 2 * T2;
T1 = T1 + T1dot * t ;
T2 = T2 + T2dot * t ;
% Dynamic equations
xddot = ( T1 + T2 ) * sin(theta) / m;
yddot = ( T1 + T2 ) * cos(theta) / m - g;
thetaddot = L*(T1 + T2)/I;
% Integrate velocity
xdot = xdot + xddot * t;
ydot = ydot + yddot * t;
thetadot = thetadot + thetaddot * t;
% % Position
% x = x + xdot * t;
% y = y + ydot * t;
% theta = theta + thetadot * t;
output = [xdot     ;...
          ydot     ;...
          thetadot ;...
          xddot    ;...
          yddot    ;...
          thetaddot;...
          T1dot    ;...
          T2dot];
end

