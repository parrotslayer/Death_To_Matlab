
%%% Base file for AMME5520 Assignment 2
% Navigation and control through an obstacle field. The objective is to go
% from left to right.
clear
clc
close all
Width = 300;
Height = 150;

% Randomly generate some obstacles (number, average size as parameters)
numObst = 20; % Control the number of obstacles, e.g. 10, 20, 30, 40.
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
num_nodes = 2000;
num_edgenodes = 10;
nodes = PRM_node_generator(As,cs,num_nodes,num_edgenodes);
no_of_neighbours = 10;
costplot = zeros(num_nodes,num_nodes);
cost = NaN(num_nodes,num_nodes);
for index1 = 1:length(nodes(1,:))
    node = nodes(:,index1);
    [idx,dist] = knnsearch(nodes',node','k',no_of_neighbours);
    for index2 = 1:length(idx)
         costplot(idx(index2),index1) = dist(index2);
         costplot(index1,idx(index2)) = dist(index2);
         cost(idx(index2),index1) = dist(index2);
         cost(index1,idx(index2)) = dist(index2);       
    end
end
plot(nodes(1,:),nodes(2,:),'g+')
gplot(costplot,nodes');
costplotcheck = costplot;
% Check the collisions for each path
for n = 1 : num_nodes
    for m = 1 : num_nodes
        if n == m
            continue
        end
        if costplot(n,m) ~= 0
            % Cycle through the
            for k = 1 : numObst
                % If badpath, returns 1 else 0
                badpath = CheckCollision(nodes(:,n),nodes(:,m),As{k},cs{k});
                if badpath
                    costplotcheck(n,m) = 0;
                    break
                end
            end 
        end
    end
end
% Plot Corrected Version with no intersecting paths
figure
for k = 1 : numObst
     Ellipse_plot(As{k},cs{k});
end
plot(nodes(1,:),nodes(2,:),'g+')
gplot(costplotcheck,nodes');
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');

djikstrasinput = costplotcheck;
for n = 1 : num_nodes
    for m = 1 : num_nodes
        if djikstrasinput(n,m) == 0
            djikstrasinput(n,m) = NaN;
        end
    end
end
[~,beginpoint] = min(nodes(1,:));
end_index = find(nodes(1,:)>250);
% [~,endpoint] = max(nodes(1,:));
% for m = 1 : length(end_cand)
[optimalpath,cost] = Djikstras(beginpoint,end_index,djikstrasinput);
[best_cost,best_path] = min(cost);
best_path_route = optimalpath{best_path};

%    optimalpath_store{m} = optimalpath;
%    cost_store(m) = cost;
% end

% minroute = find(cost_store == min(cost_store));
% final_route = optimalpath_store{minroute};
for n = 1 : length(best_path_route)
    pathnodes(:,n) = nodes(:,best_path_route(n));
end
pathbegin = [0;pathnodes(2,1)];
pathend = [300;pathnodes(2,end)];
pathnodes = [pathbegin,pathnodes,pathend];
figure
for k = 1 : numObst
     Ellipse_plot(As{k},cs{k});
end
plot(pathnodes(1,:),pathnodes(2,:),'-x');
plot([dimensions(1) dimensions(2)],[dimensions(3) dimensions(3)],'r');
plot([dimensions(1) dimensions(2)],[dimensions(4) dimensions(4)],'r');
grid on
grid minor
