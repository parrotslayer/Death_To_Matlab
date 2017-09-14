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

im = imread(image_R{1,1});

imshow(im)



