clc
clear all
close all

lego = imread('legomontage3.png');
%lego = imread('legobricks012.jpg');


figure
imshow(lego)
title('Image Loaded')

% Apply Filter
%[BW,maskedRGBImage] = Im1_Blue(lego);
[BW,maskedRGBImage] = Montage_LGreen2(lego);

figure
imshow(BW)

figure
imshow(maskedRGBImage)
title('Filtered Image')

%% Get info on regions
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);
figure
imshow(BW2)

BW3 = bwareafilt(BW,1);
figure
imshow(BW3)

%BW2 = bwareafilt(BW,1);


regions_filt = regionprops(BW3, 'Centroid', 'Area');