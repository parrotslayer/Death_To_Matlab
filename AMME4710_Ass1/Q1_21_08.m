%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     I = ro*N*s    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% where:
% ro -> albedo (reflectance) - 
% N -> Normal
% s -> direction of the light

% im_array - has all the images

%only first pixel for now


% I = s*g
% b = Ax    ...  x = A^-1 b  or  A^-1\b

close all
clear all
clc

load('facedata_yaleB01.mat')

[rows,cols,N] = size(im_array);

A = light_dirs;

albedo = zeros(rows,cols);
normal = zeros(rows,cols,3);

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

%check if (dp/dy - dq/dx)^2 is small

%calculate integral
height_map = (cumsum(p,2) + cumsum(q,1)) ;

% set 1st pixel to 0
offset = height_map(1,1);

%offset array
height_map = height_map - offset;

display_face_model(albedo, height_map)

%% Part 2 Outlier detection

I_est = zeros(rows,cols,N);
residuals = zeros(rows,cols,N);
outlier = zeros(rows,cols,N);
im_array2 = zeros(rows,cols,3,N);

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
        
        for k = 1:N
            if residuals(i,j,k) > upper || residuals(i,j,k) < lower
                im_array2(i,j,:,k) = cat(3,100,100,100);
                outlier(i,j,k) = 1;
            else
                im_array2(i,j,:,k) = reshape(cat(3,im_array(i,j,k),im_array(i,j,k),im_array(i,j,k)),[3,1]);                
            end               
        end
        
    end
end
   
%% show montage
figure
imshow(im_array2(:,:,:,1))
%%
figure
montage(im_array2)
%%
figure
montage(reshape(outlier,rows,cols,1,N))