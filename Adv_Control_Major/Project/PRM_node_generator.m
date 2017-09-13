function PRM_nodes = PRM_node_generator(As,cs,num_nodes)
num_obs = length(As(:,1));
num_valid_obs = 1;
point = NaN(2,num_nodes);
while num_valid_obs < num_nodes + 1
    n = num_valid_obs;
    hit = ones(1,num_obs);
    % Check if the current point is in an object
    while sum(hit) >= 1 % IF INSIDE ANY OBJECT
        %Generate point
        point(:,n) = [40 + rand * 220 ; rand * 150];
        % Check if it hits any objects
        for k = 1 : num_obs
            hit(k) = CheckCollision(point(:,n),point(:,n),As{k},cs{k});
        end
    end
    num_valid_obs = num_valid_obs + 1;
end
PRM_nodes = point;

