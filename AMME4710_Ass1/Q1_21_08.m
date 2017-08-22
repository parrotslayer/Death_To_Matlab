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

load('facedata_yaleB07.mat')

[h,w,N] = size(im_array);

A = light_dirs;

albedo = zeros(h,w);
normal = zeros(h,w,3);

for i = 1:h   
    for j = 1:w
     
    b = reshape(im_array(i,j,:),[N,1]);     % 1x1x64 -> 64x1
    x = A\(double(b)/255);                  % Intensity b (uint8) scaled to [0-1]
    
    %compute albedo and normals
    albedo(i,j) = norm(x);
    nhat = x/albedo(i,j);
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
%height_map = height_map - offset;

display_face_model(albedo, height_map)

%figure
%imagesc(albedo)
%colormap(gray)