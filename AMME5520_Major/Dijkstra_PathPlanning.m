%% AMME5520 Advanced Systems
% Assignment 1 Q1
clc
clear all
close all

%*************** Starting and Ending Nodes for user to edit ***************
start = 1;      %start
ending = 2;     %end
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

D = [NaN,NaN,67.5267445282944,257.604928108002,205.510412641924,121.115889764507,NaN,NaN,NaN,NaN;NaN,NaN,184.794152867616,102.523904974127,87.3620946084542,35.1170417271381,56.5464272842381,81.6548992792183,64.5180265101702,77.2909860528746;NaN,NaN,NaN,NaN,252.939219211200,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,111.427643628032,136.604049398712,114.407507148725,81.8935858781909,NaN,NaN;NaN,NaN,NaN,NaN,NaN,111.208133335023,36.2319012936076,29.9826486137487,82.3369643027440,NaN;NaN,NaN,NaN,NaN,NaN,NaN,75.7497882810125,NaN,56.6927589813483,54.1291244607479;NaN,NaN,NaN,NaN,NaN,NaN,NaN,48.5933097464955,49.7772271559952,79.7234552958686;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,31.1491841568940;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
[size_D,cols] = size(D);
Q = zeros(1,cols);
%Q = [0,0,0,0,0,0,0];    %Open Set to be visited
Q(start) = 1;       %set starting location to 1

V = ones(1,cols)*inf;    %Min cost to each node
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
        if isnan(D(x,j)) ~= 1 
            %if there is a path, store the node into R
            paths = paths + 1;
            R(paths) = j;
        end
    end
    
    %% For each reachable node, calculate the distance from x
    [rows,cols] = size(R);
    if norm(R) > 0
        for k=1:cols
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
    
end   %end while loop

%% Now get shortset path from Start to Ending
found = 0;  %check for if solution found
x = ending;
i = 0;
inc = 0

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
    inc = inc + 1
    if inc > 1000
        break;
    end
end

%rotate path matrix
path = rot90(path,2);
disp(['Shortest Path = ' num2str(path)]);
disp(['Minimum Cost = ' num2str(V(ending))]);