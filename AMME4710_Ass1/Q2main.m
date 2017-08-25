clc
clear all
close all

lego = imread('legomontage.png');

figure
imshow(lego)
title('Image Loaded')

% Apply Filter
%[BW,maskedRGBImage] = Im1_Blue(lego);
[BW,maskedRGBImage] = Montage_Red(lego);

figure
imshow(BW)

figure
imshow(maskedRGBImage)
title('Filtered Image')
