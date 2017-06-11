function y = meas(x, variance)
%% Define measurements

y = normrnd(x,variance); %start with perfect state feedback.
