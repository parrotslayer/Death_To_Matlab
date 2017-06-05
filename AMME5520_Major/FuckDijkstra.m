%% AMME5520 Advanced Systems
% Assignment 1 Q1
clc
clear all
%close all

%*************** Starting and Ending Nodes for user to edit ***************
start = 1;      %start
ending = 2;     %end
%**************************************************************************

inf = 10e10;
inc = 0;

% distances between each node
D = [NaN,NaN,NaN,16.3497783534859,91.9115269563890,76.1367615426584;NaN,NaN,101.677699094411,158.816177309504,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,71.8839784656566;NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN];
D = [NaN,NaN,63.4426592066462,54.3109644262895,NaN,77.0401863471982,88.6844123329914;NaN,NaN,147.115339428364,116.515278193398,29.9493750371342,NaN,NaN;NaN,NaN,NaN,95.2652901879952,122.158110527274,45.6567676915738,46.9062051221530;NaN,NaN,NaN,NaN,87.1912700273064,125.145202997126,NaN;NaN,NaN,NaN,NaN,NaN,NaN,83.6305064449840;NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN];
D = [NaN,NaN,118.934247403411,92.6237800205256,NaN,106.954620713786,NaN;NaN,NaN,NaN,NaN,106.584435372121,69.0583295763258,40.4897754885081;118.934247403411,NaN,NaN,57.1112185406558,NaN,NaN,NaN;92.6237800205256,NaN,57.1112185406558,NaN,NaN,NaN,NaN;NaN,106.584435372121,NaN,NaN,NaN,173.555762390938,120.865347447525;106.954620713786,69.0583295763258,NaN,NaN,173.555762390938,NaN,60.9131199638206;NaN,40.4897754885081,NaN,NaN,120.865347447525,60.9131199638206,NaN];
[size_D,cols] = size(D);

% Initialise Q = Open Set, nodes to be visited
Q = ones(1,size_D);

% Initaialise S visited set
S = zeros(1,size_D);

% Initilise cost matrix V, store min cost to each node
V = ones(1,size_D)*inf;    %Min cost to each node
V(start) = 0;   %cost to starting node is 0

% Initialise Pi = source node
Pi = zeros(1,size_D);

% loop as long as there are nodes to reach
while norm(Q) ~= 0
    
    % Select first X from Q
    for j=1:size_D
        if Q(j) ~= 0
            x = j;
            break
        end
    end
    % Move X from Q -> S
    Q(x) = 0;
    S(x) = 1;
    
    % Find number of paths from X
    paths = 0;  %number of path from x
    R = zeros(1);   %array to store path sources
    for j = 1:size_D
        %check if there is a path. Check along the rows
        if isnan(D(x,j)) ~= 1
            %if there is a path, store the node into R
            paths = paths + 1;
            R(paths) = j;
        end
    end
    disp(x)
    
    if norm(R) > 0
        [rows,num_paths] = size(R);
        % Repeat for all y Paths from X
        for k = 1:num_paths
            y = R(k);   % y is the current neighbour we are looking at
            dist = V(x) + D(x,y);  % calc distance            
            %check if distance in this path is less than current path
            if dist < V(y)
                V(y) = dist; %update distance
                Pi(y) = x;    %update Pi
                Q(y) = 1;     %put Y into open set as a new shortest path found 
            end
        end % end for loop
    end % end norm(R) > 0
end   %end while loop
