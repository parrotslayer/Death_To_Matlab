%% AMME5520 Advanced Systems
% Assignment 1 Q1
clc
clear all
close all

%*************** Starting and Ending Nodes for user to edit ***************
start = 1;      %start at A
ending = 6;     %end at F
%**************************************************************************

size_D = 7;
inf = 10e10;

% distances between each node
D = [0 10 1 0 0 0 0;
    0 0 5 0 0 0 0;
    2 5 0 3 4 0 0;
    0 0 0 0 0 0 2;
    0 0 0 0 0 2 1;
    0 3 0 0 0 0 0;
    0 0 0 2 0 5 0];

Q = [0,0,0,0,0,0,0];    %Open Set to be visited
Q(start) = 1;       %set starting location to 1

V = ones(1,7)*inf;    %Min cost to each node
V(start) = 0;   %cost to starting node is 0

% loop as long as there are nodes to reach
while norm(Q) ~= 0
    %x is the first non zero term in the Q matrix
    for j=1:size_D
        if Q(j) ~= 0
            x = j;
        end
    end
    
    %% build R, reachable nodes from x
    paths = 0;  %number of path from x
    R = zeros(1);   %array to store path sources
    for j = 1:size_D
        %check if there is a path
        if D(x,j) ~= 0
            %if there is a path, store the node into R
            paths = paths + 1;
            R(paths) = j;
        end
    end
    
    %% For each reachable node, calculate the distance from x
    [rows,cols] = size(R);
    for k=1:cols
        %get distance from x to kth node in R
        dist = V(x) + D(x,R(k));
        %check if distance in this path is less than current path
        if dist < V(R(k))
            V(R(k)) = dist; %update distance
            Q(R(k)) = 1;    %update Q
        end
    end
    
    Q(x) = 0;   %finished looking at x
    
end   %end while loop

%% Now get shortset path from Start to Ending
found = 0;  %check for if solution found
x = ending;
i = 0;

while found == 0
    for j=1:size_D
        test = V(x) - D(j,x);   %distance from x to j
        %if the distance from x to j is the shortest, then is a path
        if test == V(j) && D(j,x) ~= 0
            i = i + 1;
            path(i) = x;    %store node into path
            x = j;  %new node is j
            break
        end
    end
    %if we have reached the starting node, stop
    if x == start
        found = 1;
        path(i+1) = start;
    end    
end

%rotate path matrix
path = rot90(path,2);
disp(['Shortest Path = ' num2str(path)]);
disp(['Minimum Cost = ' num2str(V(ending))]);