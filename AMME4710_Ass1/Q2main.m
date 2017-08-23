clc
clear all
close all

lego = imread('legobricks007.jpg');

figure
imshow(lego)
title('Image Loaded')

% Apply Filter
[BW,maskedRGBImage] = Im1_Blue(lego);

%figure
%imshow(BW)

figure
imshow(maskedRGBImage)
title('Filtered Image')
