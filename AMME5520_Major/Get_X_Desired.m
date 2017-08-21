function X_Desired = Get_X_Desired(path,spacing,vel)
num_points = length(path(1,:));
direction = NaN(2,num_points-1);
output = [];
output2 = [];
spacedinterval = [];
counter = 0;
spaceddirection = [];
for n = 1 : num_points-1
    counter = counter + 1;
    X1 = path(:,n);
    X2 = path(:,n+1);
    direction(:,n) = (X2-X1)/norm(X2-X1);
    X = X1;
    while X(1) < X2(1)  
    X = X + direction(:,n) * spacing;
        if X(1) > X2(1)
            X = X2;
        end
         spacedinterval = [spacedinterval,X];
         spaceddirection = [spaceddirection,direction(:,counter)*vel];
    end    
end
output2 = [spaceddirection,[0;0]];
spacedinterval = [path(:,1),spacedinterval];
otherstates = zeros(1,length(output2(1,:)));
thrust = ones(1,length(output2(1,:))) * 0.612 * 9.81/2;
X_Desired = [spacedinterval;otherstates;output2;otherstates;thrust;thrust];
end