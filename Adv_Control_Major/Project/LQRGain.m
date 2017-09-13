function gain = LQRGain(modelparams)
m = modelparams(1);
g = modelparams(2);
L = modelparams(3);
I = modelparams(4);

%% Part D - Get the K matrix
% Defining the A Matrix 

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

% The outputs we want to observe are X, Y and theta 
C = [1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0];

% I am placing a weight of one for each of the variables that deviate from
% the equilibrium (X,Y and theta) as seen by the fact only the first three
% rows contain ones (the weights) to the relevant elements
Q = C' * C;
w1 = 0.25^2;
w2 = 0.25^2;
w3 = 0.20^2;

% Apply Bryson's rule to add weights to x, y, theta
Q(1,1) = 1/w1;
Q(2,2) = 1/w2;
Q(3,3) = 1/w3;

% Again placing a weight of one on the control effort (subject to change
% later)
R = [1 0;
    0 1];

[K,~] = lqr(A,B,Q,R);

gain = K;
end