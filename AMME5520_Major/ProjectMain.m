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
ending_point = [100,120];

% starting_point = [250,101];
% ending_point = [100,101];

for k = 1:numObst
    % Generate ellipse in the form (x-c)'A(x-c)=1
    L = randn(2,2);
    As{k} = (0.4*eye(2)+ L'*L)/Adim^2;
    tmp = rand(2,1);
    cs{k} = [leftlim+(rightlim-leftlim)*tmp(1);dimensions(3)+(dimensions(4)-dimensions(3))*tmp(2)];
    
       %  As{k} = [0.0289, 0.0036; 0.0036, 0.0052];
       %  cs{k} = [170.3786;100.7861];
    
    %Add offset to ellipse, is proportional
    Ass(:,:,k) = As{k};
    Ass(1,1,k) = Ass(1,1,k)*buffer;
    Ass(1,2,k) = Ass(1,2,k)*buffer;
    Ass(2,1,k) = Ass(2,1,k)*buffer;
    Ass(2,2,k) = Ass(2,2,k)*buffer;
           
    Ellipse_plot(As{k},cs{k},'b');
    Ellipse_plot(Ass(:,:,k),cs{k},'r');
    hold on
    
    hit = CheckCollision(starting_point,ending_point,Ass(:,:,k),cs{k});
    
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
N = 100; %number of points/nodes required
Width = 300;
Height = 150;
%Ass(:,:,k)
%cs{k}
starting_point = [250,120];
ending_point = [100,80];

%****************************** FUNCTON BEGINS **************************
% Constants and variables init
params = ComputePath(Ass,cs);

G = NaN(N,2); %pre-allocate for speed
Distances = NaN(N); % pre-allocate for speed
num_lines = 0;

[r,c,num_obs] = size(Ass);  %get number of obstacles
i = 2;  %counter starts at 2 because have start and finish point already
K = 3;  %number of nearest nodes to connect to

drawrealtime = 0;

%Init G with the start and ending points
G(1,:) = starting_point;
G(2,:) = ending_point;

% Plot Ellipses, Nodes and Paths
figure
for k = 1:numObst
Ellipse_plot(As{k},cs{k},'b');
Ellipse_plot(Ass(:,:,k),cs{k},'r');
hold on
end
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
plot(G(1,1),G(1,2),'ro')
plot(G(2,1),G(2,2),'ro')
hold on

% Begin PRM
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
        % add alpha to the last row of G
        G(i,:) = alpha;
        
        % Plot this new point
        if drawrealtime == 1
            plot(G(i,1),G(i,2),'b.')
            hold on
            title(sprintf('PRM Generation: Node %d',i))
            drawnow
        end
        %Find K nearest nodes
        %disp('Current Node = ')
        %disp(i)
        
        % Use KNNsearch to find closest distance to K nodes
        [idx,d] = knnsearch(G, G, 'k', K + 1);
        % idx returns the row number of the node
        idx = idx(:,2:K+1);
        % d is the euclidian distance arraged from shortest -> longest
        d = d(:,2:K+1);
        
        %disp(idx)
        %disp(d)
        
        %for each neighbour
        for j = 1:K
            %Check if alpha goes through an obstacle
            point1 = alpha;
            % ending point is the jth closest neighbour to alpha
            point2 = G(idx(i,j),:);
            % Check if path crosses between an obstacle
            obs_hit = 0;
            for k = 1:num_obs
                % see if there is a hit
                hit = CheckCollision(point1,point2,Ass(:,:,k),cs{k});
                % if there is a hit exit checking obstacles
                if hit == 1
                    obs_hit = 1;
                    break
                end
            end
            
            % if the path between alpha and neightbours are clear
            if obs_hit == 0
                % add the path to q/Distances if not exists
                Distances(idx(i,j),i) = d(i,j);
                % Plot the lines/path
                
                if drawrealtime == 1
                    X_realtime = [point1(1),point2(1)];
                    Y_realtime = [point1(2),point2(2)];
                    plot(X_realtime,Y_realtime,'k');
                    hold on
                    drawnow
                else
                    if isnan(point1(1)) ~= 1 && isnan(point2(1)) ~= 1
                        num_lines = num_lines + 1;
                        X_lines(num_lines,:) = [point1(1),point2(1)];
                        Y_lines(num_lines,:) = [point1(2),point2(2)];
                    end
                end
                
            end %end obs_hit == 0

        end % end for 1:K for each neightbour
        
        
    end %end if obs_hit
    
    %i = N;
end

if drawrealtime == 0
    plot(G(:,1),G(:,2),'b.')
    hold on
    for j = 1:num_lines
        plot(X_lines(j,:),Y_lines(j,:),'k')
        hold on
    end
end

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



