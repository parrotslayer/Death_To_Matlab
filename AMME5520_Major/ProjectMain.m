%%% Base file for AMME5520 Assignment 2
% Navigation and control through an obstacle field. The objective is to go
% from left to right.
clc
close all
clear all

Width = 300;
Height = 150;

% Randomly generate some obstacles (number, average size as parameters)
numObst = 3; % Control the number of obstacles, e.g. 10, 20, 30, 40.
Adim = 15; % Control the "average" size of the obstacles.

As = cell(numObst,1);
cs = cell(numObst,1);

figure(1);

dimensions = [0 Width 0 Height];

% Create obstacles only in the middle 60%
leftlim = 0.2*(dimensions(2)-dimensions(1));
rightlim = 0.8*(dimensions(2)-dimensions(1));

buffer = 0.5;    %buffer zone in percent
% Points in form [X,Y]
starting_point = [250,120];
ending_point = [100,80];

for k = 1:numObst
    % Generate ellipse in the form (x-c)'A(x-c)=1
    L = randn(2,2);
    As{k} = (0.4*eye(2)+ L'*L)/Adim^2;
    tmp = rand(2,1);
    cs{k} = [leftlim+(rightlim-leftlim)*tmp(1);dimensions(3)+(dimensions(4)-dimensions(3))*tmp(2)];
    
    %     As{k} = [0.0289, 0.0036; 0.0036, 0.0052];
    %     cs{k} = [170.3786;100.7861];
    
    %Add offset to ellipse, is proportional
    Ass(:,:,k) = As{k};
    Ass(1,1,k) = Ass(1,1,k) + Ass(1,1,k)*(buffer);
    Ass(1,2,k) = Ass(1,2,k) + Ass(1,2,k)*(buffer);
    Ass(2,1,k) = Ass(2,1,k) + Ass(2,1,k)*(buffer);
    Ass(2,2,k) = Ass(2,2,k) + Ass(2,2,k)*(buffer);
    
    Ellipse_plot(As{k},cs{k});
    Ellipse_plot(Ass(:,:,k),cs{k});
    
    hit = CheckCollision(starting_point,ending_point,Ass(:,:,k),cs{k});
    disp(hit);
    
end
X = [starting_point(1),ending_point(1)];
Y = [starting_point(2),ending_point(2)];
plot(X,Y,'k');
hold on
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
hold off

%% Students add code here for planning
% - Students to modify ComputePath, to include collision checker, PRM, and
% path planner.

%****************************** Input Variables *************************
N = 10; %number of points/nodes required
Width = 300;
Height = 150;
%Ass(:,:,k)
%cs{k}

%****************************** FUNCTON BEGINS **************************
params = ComputePath(Ass,cs);

G = zeros(N,2); %pre-allocate for speed
[r,c,num_obs] = size(Ass);  %get number of obstacles
i = 0;  %counter
K = 2;  %number of nearest nodes to connect to
while i < N
    
    % generate a random point alpha = [X, Y]
    alpha(1) = rand*Width;
    alpha(2) = rand*Height;
    
    % Check if point if within any of the obstacles
    obs_hit = 0;
    for k = 1:num_obs
        % see if there is a hit
        hit = CheckCollisionPoint(alpha,Ass(:,:,k),cs{k});
        % if there is a hit
        if hit == 1
            obs_hit = 1;
            break
        end
    end
    
    if obs_hit == 0
        %increment counter i
        i = i+1;
        % add alpha to G
        G(i,:) = alpha;
        % Try to connect to K nearest nodes
        if i > K
            %Find K nearest nodes
            disp('Find Nearest Nodes and Check')
            %Check if through obstacle
            
            %If node clear
            
            %Add q to G
            
            
        end % end i > K connecting to nearest nodes
                
    end %end if obs_hit
    
    %i = N;
end

% Plot Ellipses, Nodes and Paths
figure
Ellipse_plot(As{k},cs{k});
Ellipse_plot(Ass(:,:,k),cs{k});
hold on
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
hold on
% Node G = [X, Y]
plot(G(:,1),G(:,2),'b.')

title('Compute Path')

%***************************** END **************************************


%% Simulate Closed-loop system
% students to modify functions

h = 0.01; % 100Hz sample time. Modify as desired.

x0 = zeros(8,1);  %Initial condition.

stop = 0;
ts = 0;
xs = x0;

% This variable can be used to collect any parameters for simulation
dynparams = []


k = 1;
while (stop ~= 1)
    % Current state.
    xt = xs(:,k);
    
    % Get current measurement, compute control.
    yt = meas(xt);
    [ut, params] = ComputeControl(yt,params);
    
    
    % Use Runge Kutta to approximate the nonlinear dynamics over one time
    % step of length h.
    xs(:,k+1) = RungeKutta4(@QuadDynamics, xt, ut, 0, h, dynparams);
    
    ts(k+1) = ts(k)+h;
    us(:,k) = ut;
    
    if (k>100)
        stop = 1;  % Stop after 100 time steps. You will need to change this.
        % Should stop after goal reached or collision.
    end
    k = k+1;
end



