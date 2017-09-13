close all
clear all
clc

load('facedata_yaleB07.mat')

load(['..','assignment2_stereodata',filesep,'camera_pose_data.mat']) ;

% load(fullfile('..', 'folderName1', 'folderName2', 'filename.ext'))

[rows,cols,N] = size(im_array);

A = light_dirs;


%% Get surface normals
%preallocate for speed
albedo = zeros(rows,cols);
normal = zeros(rows,cols,3);
p = zeros(rows,cols);
q = zeros(rows,cols);

% loop for each pixel
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
