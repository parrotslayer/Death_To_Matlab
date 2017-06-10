% inputxy as a 2x44 vector of nodes
% spacing = dt*vel
% vel = velocity = 3m/s

function path = Get_X_desired(inputxy,spacing,vel, dynparams)
m = dynparams(1);
g = dynparams(2);
L = dynparams(3);
I = dynparams(4);
freq = dynparams(5);

num_points = length(inputxy(1,:));
direction = NaN(2,num_points-1);

spaced_XY = [];
counter = 0;
spaced_XY_dot = [];

% repeat for the number of points
for n = 1 : num_points-1
    counter = counter + 1;
    X1 = inputxy(:,n);
    X2 = inputxy(:,n+1);
    direction(:,n) = (X2-X1)/norm(X2-X1);
    X = X1;
    while X(1) < X2(1)  
    X = X + direction(:,n) * spacing;
        if X(1) > X2(1)
            X = X2;
        end
         spaced_XY = [spaced_XY,X];
         spaced_XY_dot = [spaced_XY_dot,direction(:,counter)*vel];
    end    
end

% Add final values to position and velocity
output2 = [spaced_XY_dot,[0;0]];
spaced_XY = [inputxy(:,1),spaced_XY];

% Create array of zeros for th, th'
zero_states = zeros(1,length(output2(1,:)));

% T1 and T2 are just mg/2
thrust = ones(1,length(output2(1,:)))*m*g/2;

% Combine arrays together
path = [spaced_XY;zero_states;output2;zero_states;thrust;thrust];
end