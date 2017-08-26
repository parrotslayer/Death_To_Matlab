clc
clear all
close all

% Load the image
%lego = imread('legomontage3.png');
lego_num = ['001';'002';'003';'004';'005';'006';'007';'008';'009';'010';...
    '011';'012';'013';'014';'015';'016';'017'];
%choose the image to display
N = 1;

for N = 1:17

filename = ['legobricks',lego_num(N,:),'.jpg'];
lego = imread(filename);

%Preallocate arrays
% Red, Dark Green, Blue, Light Green, Yellow, Orange
All_centroids = NaN(6,2);
All_boxes = NaN(6,4);

figure
imshow(lego)
title(['Image Loaded: ',filename])

%% Red Filtering
% Apply filter
[BW,Red] = Montage_Red(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take largest block (no conflicts)
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_Red = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a red block detected
[n,trump] = size(regions_Red);
if n > 0
    
    %Store information in arrays
    centroid = regions_Red.Centroid;
    All_centroids(1,:) = centroid;     
    box = regions_Red.BoundingBox;
    All_boxes(1,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Red';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
end

%% Dark Green Filtering
% Apply filter
[BW,DGreen] = Montage_DGreen2(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take largest block (no conflicts)
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_DGreen = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a dark green block detected
[n,trump] = size(regions_DGreen);
if n > 0
    
    %Store information in arrays
    centroid = regions_DGreen.Centroid;
    All_centroids(2,:) = centroid;     
    box = regions_DGreen.BoundingBox;
    All_boxes(2,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Dark Green';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
end

%% Blue Filtering
% Apply filter
[BW,Blue] = Montage_Blue(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take largest block, high confidence only wanted block is there
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_Blue = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a blue block detected
[n,trump] = size(regions_Blue);
if n > 0
    
    %Store information in arrays
    centroid = regions_Blue.Centroid;
    All_centroids(3,:) = centroid;     
    box = regions_Blue.BoundingBox;
    All_boxes(3,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Blue';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
end


%% Light Green Filtering
% Apply filter
[BW,LGreen] = Montage_LGreen2(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take largest block (no conflicts)
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_LGreen = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a light green block detected
[n,trump] = size(regions_LGreen);
if n > 0
    
    %Store information in arrays
    centroid = regions_LGreen.Centroid;
    All_centroids(4,:) = centroid;     
    box = regions_LGreen.BoundingBox;
    All_boxes(4,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Light Green';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
end
%% Yellow Filtering
% Apply filter
[BW,Yellow] = Montage_Yellow1(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take several largest blocks cos can be lgreen, yellow or orange
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_Yellow = regionprops(BW3, 'Centroid', 'Area','BoundingBox');

%check if there is a yellow block detected
[n,trump] = size(regions_Yellow);
if n > 0
    
    %Store information in arrays
    centroid = regions_Yellow.Centroid;
    All_centroids(4,:) = centroid;     
    box = regions_Yellow.BoundingBox;
    All_boxes(4,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Yellow';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
end

end