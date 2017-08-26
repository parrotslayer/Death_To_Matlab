clc
clear all
close all

% Load the image
%lego = imread('legomontage3.png');
lego = imread('legobricks009.jpg');

figure
imshow(lego)
title('Image Loaded')

box_length = 200;
centroid = [1000,2000];
% rectangle [x y w h] 
rectangle('Position',[centroid(1)-box_length/2,centroid(1)-box_length/2,box_length,box_length],'LineWidth',6,'EdgeColor','r');
string = 'color';
text(centroid(1),centroid(2),string)
%% Blue Filtering
% Apply filter
[BW,Blue] = Montage_Blue(lego);

figure
imshow(Blue)
title('Post Filter Blue')

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

figure
imshow(BW2)
title('Post Thresholding BW')

% Take largest block (no conflicts)
BW3 = bwareafilt(BW,1);
figure
imshow(BW3)
title('Show largest Region')

% Details on the final block
regions_blue = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%Create bounding box
box = regions_blue.BoundingBox;
rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
string = 'Blue';
text(box(1),box(2)+100,string,'Color','White','FontSize',14)
%% Dark Green Filtering
% Apply filter
[BW,DGreen] = Montage_DGreen2(lego);

figure
imshow(DGreen)
title('Post Filter Dark Green')

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

figure
imshow(BW2)
title('Post Thresholding BW')

% Take largest block (no conflicts)
BW3 = bwareafilt(BW,1);
figure                                                                                                                                                                                                                                                                                                                                          
imshow(BW3)
title('Show largest Region')

% Details on the final block
regions_DGreen = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%Create bounding box
box = regions_DGreen.BoundingBox;
rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
string = 'Dark Green';
text(box(1),box(2)+100,string,'Color','White','FontSize',14)
%% Red Filtering
% Apply filter
[BW,Red] = Montage_Red(lego);

figure
imshow(Red)
title('Post Filter Red')

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

figure
imshow(BW2)
title('Post Thresholding BW')

% Take largest block (no conflicts)
BW3 = bwareafilt(BW,1);
figure
imshow(BW3)
title('Show largest Region')

% Details on the final block
regions_Red = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%Create bounding box
box = regions_Red.BoundingBox;
rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
string = 'Red';
text(box(1),box(2)+100,string,'Color','White','FontSize',14)
%% Light Green Filtering
% Apply filter
[BW,LGreen] = Montage_LGreen2(lego);

figure
imshow(LGreen)
title('Post Filter Light Green')

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

figure
imshow(BW2)
title('Post Thresholding BW')

% Take largest block (no conflicts)
BW3 = bwareafilt(BW,1);
figure
imshow(BW3)
title('Show largest Region')

% Details on the final block
regions_Lgreen = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%Create bounding box
box = regions_Lgreen.BoundingBox;
rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
string = 'Light Green';
text(box(1),box(2)+100,string,'Color','White','FontSize',14)
%% Yellow Filtering
% Apply filter
[BW,Yellow] = Montage_Yellow1(lego);

figure
imshow(Yellow)
title('Post Filter Yellow')

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

figure
imshow(BW2)
title('Post Thresholding BW')

% Take largest block (no conflicts)
BW3 = bwareafilt(BW,1);
figure
imshow(BW3)
title('Show largest Region')

% Details on the final block
regions_Yellow = regionprops(BW3, 'Centroid', 'Area','BoundingBox');

%Create bounding box
box = regions_Yellow.BoundingBox;
rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
string = 'Yellow';
text(box(1),box(2)+100,string,'Color','White','FontSize',14)

