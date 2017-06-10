%% AMME5520 Advanced Systems
function [min_dist, shortest_path] = Dijkstra_Path(D,start,ending)
%% Init Variables

% Get size of data
[size_D,cols] = size(D);

inf = 10e10;

% Initialise Q = Open Set, nodes to be visited
Q = ones(1,size_D);

% Initaialise S visited set
S = zeros(1,size_D);

% Initilise cost matrix V, store min cost to each node
V = ones(1,size_D)*inf;    %Min cost to each node
V(start) = 0;   %cost to starting node is 0

% Initialise Pi = source node
Pi = zeros(1,size_D);

%% Find Shortest Distances to Each Path
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

%% Find Shortest Path
i = 1;
steps = 0;
found = 0;
Path = ending;
x = ending;
while found == 0
    try
        i = i + 1;
        steps = steps + 1;
        Path(i) = Pi(x);
        x = Pi(x);
        if x == start
            found = 1;
        end
    % if the obstacles are in bad places no possible solution    
    catch
        disp('WARNING: Could not find path solution!')
        disp('Check start and destination positions do not lie within an obstacle')
        disp('Also check if obstacles dont block off a path from start to finish')
        disp('Try run ProjectMain again to re-generate obstacles')
        error('Daijkstra Path Finding Has Failed!')
        
    end
end
Path = rot90(Path,2);
disp(['Shortest Path = ' num2str(Path)]);
disp(['Minimum Cost = ' num2str(V(ending))]);
min_dist = V(ending);
shortest_path = Path;

