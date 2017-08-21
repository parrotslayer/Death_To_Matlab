function [Storepaths,cost_store] = dijkstra(start,destination,map)
% initialise variable for counting the nodes
n = size(map,1);
check(1:n) = 0; % Visited points
dist(1:n) = inf; %distance to each node, from the start
prev(1:n) = NaN;
dist(start) = 0;
counter = 0;
% While the set of nodes visited is not complete
while sum(check) ~= n
    candidate = [];
    for i = 1 : n
        if check(i) == 0 % Check unvisited points
            % Incorporate the distances, into a set
           	candidate = [candidate, dist(i)];
        else % If it's already been visited, set it as infinity
            candidate = [candidate, inf];
        end
    end
    % Check for the minimum (checked already has infinite value)
    [~, u] = min(candidate) ;
    %
    check(u) = 1 ; % Mark this node as visited
    for i = 1 : n
        % Check if current distance is less than previous acquired distance
        if dist(u) + map(u,i) < dist(i) % if it is, replace it
            dist(i) = dist(u) + map(u,i); 
            prev(i) = u ; % remember the shortest path before
        end
    end
    counter = counter + 1;
    if counter > length(check)
        break
    end
end
% Go through all the destinations, find the shortest path for each one
%Initialise shortest path array, with starting position from the
%destination
for i = 1 : length(destination)
shortestp = destination(i); 
% While the first value of the shortestp matrix is not the start
while shortestp(1) ~= start
    % Go through the v, adding the previous node in the shortest 
    % path to the current node to the path.
        shortestp=[prev(shortestp(1)), shortestp];
end
Storepaths{i} =  shortestp;
cost_store(i) = dist(destination(i));
clear shortestp
end
end

