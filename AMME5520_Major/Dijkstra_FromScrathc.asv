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
D = [NaN,NaN,93.7314679622592,79.2923263410607,4.94647097693440,195.109585884427,45.8092764922902;NaN,NaN,NaN,NaN,160.180834932129,NaN,NaN;NaN,NaN,NaN,21.1581249329101,94.8757055598109,NaN,48.8593101341331;NaN,NaN,NaN,NaN,79.5855920307577,NaN,38.7509954794143;NaN,NaN,NaN,NaN,NaN,199.775165753784,47.7468515812093;NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN];

[size_D,cols] = size(D);

% Initialise Q = Open Set, nodes to be visited
Q = zeros(1,cols);
Q(start) = 1;       %set starting location to 1
% fudges
Q(ending) = 1; 


% Initilise cost matrix V, store min cost to each node
V = ones(1,cols)*inf;    %Min cost to each node
V(start) = 0;   %cost to starting node is 0

% loop as long as there are nodes to reach
while norm(Q) ~= 0
    %x is the first non zero term in the Q matrix
    for j=1:size_D
        if Q(j) ~= 0
            x = j;
            break;
            if x == 2
                disp(x)
            end
        end
    end
    
    %% build R, reachable nodes from x
    paths = 0;  %number of path from x
    R = zeros(1);   %array to store path sources
    for j = 1:size_D
        %check if there is a path
        if isnan(D(x,j)) ~= 1 
            %if there is a path, store the node into R
            paths = paths + 1;
            R(paths) = j;
        end
    end
    
    %% For each reachable node, calculate the distance from x
    [rows,num_paths] = size(R);
    % make sure there is a path to check
    if norm(R) > 0
        for k=1:num_paths
            %get distance from x to kth node in R
            try
            dist = V(x) + D(x,R(k));
            catch
                break
            end
            %check if distance in this path is less than current path
            if dist < V(R(k))
                V(R(k)) = dist; %update distance
                Q(R(k)) = 1;    %update Q
            end
        end
    end
    Q(x) = 0;   %finished looking at x
    
    inc = inc + 1;
    if inc > 1000
        break
    end
end   %end while loop
