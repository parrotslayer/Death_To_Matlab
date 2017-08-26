clc
clear all
close all

% Load the image
lego = imread('legomontage3.png');
%lego = imread('legobricks012.jpg');

figure
imshow(lego)
title('Image Loaded')

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
regions_blue = regionprops(BW3, 'Centroid', 'Area');

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
regions_DGreen = regionprops(BW3, 'Centroid', 'Area');

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
regions_Red = regionprops(BW3, 'Centroid', 'Area');

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
regions_Lgreen = regionprops(BW3, 'Centroid', 'Area');