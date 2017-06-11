function y = meas(x, variance)
%% Define measurements


%y = x; %start with perfect state feedback.

y = normrnd(x,variance);
