%Q1
clear all
close all
clc

load('terrain.mat')
load('stereo_calib.mat')
load('camera_pose_data.mat')

%directory of where the images be w.r.t folder path
images_left_dir = (['assignment2_stereodata',filesep,'images_left',filesep]) ;
images_right_dir = (['assignment2_stereodata',filesep,'images_right',filesep]) ;


%% Loop
i = 1;
num_images = 49;
%World_Points = NaN(600,3,num_images);
%colors = NaN(600,3,num_images);
for i = 1:49
    clear Cam1_Points
    clear filename
    
    disp(['Current Image is Number: ',num2str(i)])
    
    image_L = strcat(images_left_dir, camera_poses.left_images(i));
    image_R = strcat(images_right_dir, camera_poses.right_images(i));
    
    %images stored in a cell for some reason, want arrays
    im1_color = imread(image_L{1,1});
    im2 = imread(image_R{1,1});
    
    % convert to greyscale
    im1 = rgb2gray(im1_color);
    % the right images are in greyscale for some reason?
    %im2 = rgb2gray(im2);
    
    %% T5Q2 Using SURF to match features
    
    % compute matched set of SURF across 2 images
    points1 = detectSURFFeatures(im1,'MetricThreshold',500);
    points2 = detectSURFFeatures(im2,'MetricThreshold',500);
    
    % Extract features out 
    [descriptors1, points1] = extractFeatures(im1, points1);
    [descriptors2, points2] = extractFeatures(im2, points2);
    
    % Matching features between two images 
    % MaxRatio default is 0.6
    [matched_pairs,trump] = matchFeatures(descriptors1, descriptors2,'MaxRatio',0.1);
    points1_matched = points1(matched_pairs(:, 1), :);
    points2_matched = points2(matched_pairs(:, 2), :);
    
    %plot correspondences    
%     figure
%     showMatchedFeatures(im1,im2,points1_matched,points2_matched,'Parent',axes);
%     title(['Correspondences for Image: ',num2str(i)])
    
    % first undistort images
    %matchedPoints1 = undistortPoints(points1_matched.Location,stereoParams.CameraParameters1);
    %matchedPoints2 = undistortPoints(points2_matched.Location,stereoParams.CameraParameters2);
    
    % dont bother undistorting images because takes too long
    matchedPoints1 = points1_matched;
    matchedPoints2 = points2_matched;
    
    % compute inlier correspondences
    % dist threshold is 0.1 by default
    [F,inliersIndex] = estimateFundamentalMatrix(matchedPoints1,...
        matchedPoints2,'Method','MSAC','NumTrials',2000, 'DistanceThreshold', 0.02);
        
    %triangulate to get point cloud wrt camera 1
    Cam1_Points = triangulate(matchedPoints1(inliersIndex,:),matchedPoints2(inliersIndex,:),stereoParams);
    
    % Convert to world coordinates for each point
    [r,c] = size(Cam1_Points);
    for j = 1:r
        World_Points(j,:,i) = camera_poses.R(:,:,i) \ (Cam1_Points(j,:)' - camera_poses.t(:,i));
        
        inliers1 = points1_matched.Location(inliersIndex,:);
        x(j,:) = round(inliers1(j,1)); y(j,:) = round(inliers1(j,2));
        
        colors(j,:,i) = (double(im1_color(y(j,:),x(j,:),:))/255);
    end
end
%%  plot camera
% figure
% plot3(Cam1_Points(:,1),Cam1_Points(:,3),-Cam1_Points(:,2),'.')
% hold on
% cam1 = plotCamera('Location',[0 0 0], 'Orientation', [1,0,0; 0,0,-1; 0,1,0],'size', 0.05);
% axis equal
% xlabel('x')
% ylabel('z')
% zlabel('y')

%% Convert from camera 1 to world coordinates
figure
for i = 1:49
    scatter3(World_Points(:,1,i), World_Points(:,2,i), -World_Points(:,3,i),3,colors(:,:,i),'.')
    hold on
end
axis equal
title('World Coordinates')
xlabel('X')
ylabel('Y')
zlabel('Z')
ylim([1000 1009])
xlim([994 1000])
zlim([-4.6, -3.5])

%% plot mesh grid
% figure
% surf(X,Y,height_grid)
% 
% title('World Coordinates')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
