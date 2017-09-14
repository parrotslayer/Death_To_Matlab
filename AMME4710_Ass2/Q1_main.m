%Q1 
clear all
close all
clc

load('terrain.mat')
load('stereo_calib.mat')
load('camera_pose_data.mat')


images_left_dir = (['assignment2_stereodata',filesep,'images_left',filesep]) ;
images_right_dir = (['assignment2_stereodata',filesep,'images_right',filesep]) ;

clear filename
image_L = strcat(images_left_dir, camera_poses.left_images(1));
image_R = strcat(images_right_dir, camera_poses.right_images(1));

im1 = imread(image_L{1,1});
im2 = imread(image_R{1,1});

% convert to greyscale
im1 = rgb2gray(im1);
% the right images are in greyscale for some reason?
%im2 = rgb2gray(im2);

%% T5Q2 Copy Pasted

% compute matched set of SURF across 2 images
points1 = detectSURFFeatures(im1,'MetricThreshold',1000);
points2 = detectSURFFeatures(im2,'MetricThreshold',1000);

[descriptors1, points1] = extractFeatures(im1, points1);
[descriptors2, points2] = extractFeatures(im2, points2);
% MaxRatio default is 0.6
[matched_pairs,trump] = matchFeatures(descriptors1, descriptors2,'MaxRatio',0.4);
points1_matched = points1(matched_pairs(:, 1), :);
points2_matched = points2(matched_pairs(:, 2), :);

%plot correspondences
%showMatchedFeatures(im1,im2,points1_matched,points2_matched)

figure
showMatchedFeatures(im1,im2,points1_matched,points2_matched,'Parent',axes);

% first undistort images
%matchedPoints1 = undistortPoints(points1_matched.Location,stereoParams.CameraParameters1);
%matchedPoints2 = undistortPoints(points2_matched.Location,stereoParams.CameraParameters2);

% dont bother undistorting images because takes too long
matchedPoints1 = points1_matched;
matchedPoints2 = points2_matched;

% compute inlier correspondences
% dist threshold is 0.1 by default
fRANSAC = estimateFundamentalMatrix(matchedPoints1,...
    matchedPoints2,'Method','MSAC','NumTrials',5000, 'DistanceThreshold', 0.01);

worldPoints = triangulate(matchedPoints1,matchedPoints2,stereoParams);

%%  plot 
figure
plot3(worldPoints(:,1),worldPoints(:,3),-worldPoints(:,2),'.')
hold on
%cam1 = plotCamera('Location',[0 0 0], 'Orientation', [1,0,0; 0,0,-1; 0,1,0],'size', 200);
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
% ylim([-1000 6000])
% xlim([-4000 4000])
% zlim([-1000, 1000])
