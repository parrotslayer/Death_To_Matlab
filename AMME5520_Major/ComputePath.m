function [min_dist, shortest_path_coordinates] = ComputePath(drawrealtime,N,K,Width,Height,dimensions,As,Ass,cs,starting_point,ending_point)
% Constants and variables init

G = NaN(N,2); %pre-allocate for speed
Distances = NaN(N); % pre-allocate for speed
num_lines = 0;

[r,c,num_obs] = size(Ass);  %get number of obstacles
i = 2;  %counter starts at 2 because have start and finish point already

%Init G with the start and ending points
G(1,:) = starting_point;
G(2,:) = ending_point;

% Plot Ellipses, Nodes and Paths
figure
for k = 1:num_obs
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
        
        % Use KNNsearch to find closest distance to K nodes
        [idx,d] = knnsearch(G, G, 'k', K + 1);
        % idx returns the row number of the node
        idx = idx(:,2:K+1);
        % d is the euclidian distance arraged from shortest -> longest
        d = d(:,2:K+1);
        
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
                Distances(i,idx(i,j)) = d(i,j);
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
    
end

if drawrealtime == 0
    plot(G(:,1),G(:,2),'b.')
    hold on
    for j = 1:num_lines
        plot(X_lines(j,:),Y_lines(j,:),'k')
        hold on
    end
end

%% Run Dijkstra
[min_dist, shortest_path] = Dijkstra_Path(Distances,1,2);
% Plot shortest path in red
[r,num_paths] = size(shortest_path);
shortest_path_coordinates = zeros(num_paths,2);
shortest_path_coordinates(1,:) = starting_point;

for j = 1:num_paths-1
    point1 = G(shortest_path(j),:);
    point2 = G(shortest_path(j+1),:);
    % Store these points into a variable
    shortest_path_coordinates(j+1,:) = point2;
    X_realtime = [point1(1),point2(1)];
    Y_realtime = [point1(2),point2(2)];
    plot(X_realtime,Y_realtime,'r','linewidth',2);
    plot(point1(1),point1(2),'ro')
    plot(point2(1),point2(2),'ro')
    hold on
end

end