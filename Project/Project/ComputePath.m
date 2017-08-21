function [pathnodes,best_cost] = ComputePath(As,cs,num_nodes)
nodes = PRM_node_generator(As,cs,num_nodes);
no_of_neighbours = 10;
costplot = zeros(num_nodes,num_nodes);
cost = NaN(num_nodes,num_nodes);
numObst = length(As);
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
[optimalpath,cost] = Djikstras(beginpoint,end_index,djikstrasinput);
[best_cost,best_path] = min(cost);
best_path_route = optimalpath{best_path};
pathnodes = NaN(2,length(best_path_route));
for n = 1 : length(best_path_route)
    pathnodes(:,n) = nodes(:,best_path_route(n));
end
pathbegin = [0;pathnodes(2,1)];
pathend = [300;pathnodes(2,end)];
pathnodes = [pathbegin,pathnodes,pathend];

