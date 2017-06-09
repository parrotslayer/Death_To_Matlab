% Discretises an array of lines into smaller nodes that correspond to
% movement at a constant velocity.

%**************************** testing ******************************
% clc
% clear all
% velocity = 10;
% delta_T = 1/100;
% shortest_path_coordinates = [0,0,10,10;
%     10,10,-20,40;
%     -20,40,10,60];

function desired_XY = discretise_path(velocity, delta_T, shortest_path_coordinates)
% get size of how small to break up the path
dS = velocity*delta_T;

% get number of lines to split up
[num_lines,c] = size(shortest_path_coordinates);

% initialise the XY vector of destination to starting point
desired_XY(1,:) = [shortest_path_coordinates(1,1), shortest_path_coordinates(1,2)];

i = 1;
for k = 1:num_lines
    % determine whether to increase or decrease values in X and Y
    
    dX = shortest_path_coordinates(k,3) - shortest_path_coordinates(k,1);
    dY = shortest_path_coordinates(k,4) - shortest_path_coordinates(k,2);
    % Find angle to travel in
    theta = atan2(dY,dX);
    
    line_end = 0;
    while line_end ~= 1
        % Now divide up the current segment
        X_next = desired_XY(i,1) + dS*cos(theta);
        Y_next = desired_XY(i,2) + dS*sin(theta);
        
        % add to the desired XY matrix
        i = i+1;
        desired_XY(i,:) = [X_next, Y_next];
        
        % See if end of path is reached
        % if dX is positive, then stop when goes past maximum
        if dX > 0 && X_next > shortest_path_coordinates(k,3)
            % move onto the next line
            line_end = 1;
            % see if dX is negative
        elseif dX < 0 && X_next < shortest_path_coordinates(k,3)
            line_end = 1;
        end
    end     %while line_end ~= 1
    
end     %end for num_lines

end