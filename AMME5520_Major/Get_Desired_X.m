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

function X_desired = Get_Desired_X(velocity, delta_T, shortest_path_coordinates, dynparams)
m = dynparams(1);
g = dynparams(2);

thrust = m*g/2;

% get size of how small to break up the path
dS = velocity*delta_T;

% get number of lines to split up
[num_lines,c] = size(shortest_path_coordinates);

% initialise the X vector of destination to starting point
X_desired(1,1:4) = [shortest_path_coordinates(1,1), shortest_path_coordinates(1,2),0,0];

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
        X_next = X_desired(i,1) + dS*cos(theta);
        Y_next = X_desired(i,2) + dS*sin(theta);
        Xd_next = velocity*cos(theta);
        Yd_next = velocity*sin(theta);

        % add to the desired XY matrix
        i = i+1;
        X_desired(i,1:4) = [X_next, Y_next, Xd_next, Yd_next];
        
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

% Add last values of X_desired
X_desired(end,:) = [shortest_path_coordinates(end,3),shortest_path_coordinates(end,4)...
    0, 0];

% Set values of thrust
X_desired(:,7) = thrust;
X_desired(:,8) = thrust;

end