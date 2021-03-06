close all
clear all
clc

load('facedata_yaleB07.mat')

[rows,cols,N] = size(im_array);

A = light_dirs;

albedo = zeros(rows,cols);
normal = zeros(rows,cols,3);
p = zeros(rows,cols);
q = zeros(rows,cols);

for i = 1:rows   
    for j = 1:cols
     
    I = reshape(im_array(i,j,:),[N,1]);     % 1x1x64 -> 64x1. Pixel values from each image
    G = A\(double(I)/255);                  % Intensity b (uint8) scaled to [0-1]
    
    %compute albedo and normals
    albedo(i,j) = norm(G);
    nhat = G/albedo(i,j);
    normal(i,j,1) = nhat(1);    %x
    normal(i,j,2) = nhat(2);    %y
    normal(i,j,3) = nhat(3);    %z
    
    %compute p and q
    p(i,j) = normal(i,j,1)/normal(i,j,3);
    q(i,j) = normal(i,j,2)/normal(i,j,3);
       
    end
end

%% Integrate Average 
%calculate integral
height_map = (cumsum(p,2) + cumsum(q,1))/2 ;

% set 1st pixel to 0
offset = height_map(1,1);

%offset array
height_map = height_map - offset;

% show the 3D face
display_face_model(albedo, height_map)
title('3D reconstution of a face')

%% Across
% set first row of p from q
p2 = p;

% set first row to integral of q 
p2(:,1) = cumsum(q(:,1),2);

%calculate integral
height_map2 = cumsum(p2,2);

% set 1st pixel to 0
offset = height_map2(1,1);

%offset array
height_map2 = height_map2 - offset;

% show the 3D face
display_face_model(albedo, height_map2)
title('Integrate p across')

%% Down
% set first row of p from q
q2 = q;

% set first col to integral of p
q2(1,:) = cumsum(p(1,:),1);

%calculate integral
height_map3 = cumsum(q2,1);

% set 1st pixel to 0
offset = height_map3(1,1);

%offset array
height_map3 = height_map3 - offset;

% show the 3D face
display_face_model(albedo, height_map3)
title('Integrate q down')

%% show faces before
figure
montage(reshape(im_array,rows,cols,1,N))
title('im array before filtering')

%% Part 2 Outlier detection

I_est = zeros(rows,cols,N);
residuals = zeros(rows,cols,N);
outlier = ones(rows,cols,N);
im_array2 = zeros(rows,cols,3,N);
im_array3 = im_array;
p_filt = zeros(rows,cols);
q_filt = zeros(rows,cols);
bad_values = 0;

for i = 1:rows
    for j = 1:cols
        for k = 1:N
            %estimated light intensity
            I_est(i,j,k) = albedo(i,j)*reshape(normal(i,j,:),3,1)'*light_dirs(k,:)';
        end
        %calc residuals
        residuals(i,j,:) = double(reshape(im_array(i,j,:),[N,1])) - reshape(I_est(i,j,:),[N,1]);
        %calc std deviation and mean for thresholds
        upper = mean(residuals(i,j,:)) + 2*std(residuals(i,j,:));
        lower = mean(residuals(i,j,:)) - 2*std(residuals(i,j,:));
        
        %find residuals
        for k = 1:N
            if residuals(i,j,k) > upper || residuals(i,j,k) < lower
                %New array holding outliers
                outlier(i,j,k) = 0;
                %Set outliers in original data to be white 255
                im_array(i,j,k) = 255;
                % images in RGB with outliers highlighted in red
                im_array2(i,j,:,k) = cat(3,255,0,0);
                %remove crappy values from the copied original data
                im_array3(i,j,k) = NaN;
                % set as rediduals
                %im_array3(i,j,k) = mean(residuals(i,j,:));
                %increment number of bad values
                bad_values = bad_values + 1;
            else
                im_array2(i,j,:,k) = reshape(cat(3,im_array(i,j,k),im_array(i,j,k),im_array(i,j,k)),[3,1]);                
            end               
        end
        
        %calc new normals without residuals
        I = reshape(im_array3(i,j,:),[N,1]);     % 1x1x64 -> 64x1. Pixel values from each image
        G = A\(double(I)/255);                  % Intensity b (uint8) scaled to [0-1]
        
        %compute albedo and normals
        albedo(i,j) = norm(G);
        nhat = G/albedo(i,j);
        normal(i,j,1) = nhat(1);    %x
        normal(i,j,2) = nhat(2);    %y
        normal(i,j,3) = nhat(3);    %z
        
        %compute p and q
        p_filt(i,j) = normal(i,j,1)/normal(i,j,3);
        q_filt(i,j) = normal(i,j,2)/normal(i,j,3);
        
    end
end
   
disp(['Number of outliers detected = ' num2str(bad_values)]);
%% Integrate filtered data w outlier removal
%check if (dp/dy - dq/dx)^2 is small

%calculate integral
height_map2 = (cumsum(p_filt,2) + cumsum(q_filt,1))/2 ;

% set 1st pixel to 0
offset = height_map2(1,1);

%offset array
height_map2 = height_map2 - offset;

% show the 3D face
display_face_model(albedo, height_map2)
title('3D reconstution of a face with outliers removed')

%% show faces after
figure
montage(reshape(im_array,rows,cols,1,N))
title('im array outliers in white')

%% show im_array2
figure
montage(im_array2)
title('Faces with outliers highlighted RED')
%%
figure
montage(reshape(outlier,rows,cols,1,N))
title('show outliers in BW')