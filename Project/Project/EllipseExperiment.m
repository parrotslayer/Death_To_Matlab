%%% Base file for AMME5520 Assignment 2
% Navigation and control through an obstacle field. The objective is to go
% from left to right.
close all
clear
clc
Width = 300;
Height = 150;

% Randomly generate some obstacles (number, average size as parameters)
numObst = 1; % Control the number of obstacles, e.g. 10, 20, 30, 40.
Adim = 15; % Control the "average" size of the obstacles.


As = cell(numObst,1);
cs = cell(numObst,1);

figure(1);

dimensions = [0 Width 0 Height];

% Create obstacles only in the middle 60%
leftlim = 0.2*(dimensions(2)-dimensions(1));
rightlim = 0.8*(dimensions(2)-dimensions(1));

for k = 1:numObst
    % Generate ellipse in the form (x-c)'A(x-c)=1
    
     L = randn(2,2);
     As{k} = (0.4*eye(2)+ L'*L)/Adim^2;
     tmp = rand(2,1);
     cs{k} = [leftlim+(rightlim-leftlim)*tmp(1);dimensions(3)+(dimensions(4)-dimensions(3))*tmp(2)];
     Ellipse_plot(As{k},cs{k});
     
end

plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
hold off