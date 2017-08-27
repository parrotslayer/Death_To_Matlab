clc
clear all 
close all

I = imread('legobricks001.jpg');
imshow(I)

background = imopen(I,strel('disk',15));

% Display the Background Approximation as a Surface
figure
surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
ax = gca;
ax.YDir = 'reverse';

figure
I2 = I - background;
imshow(I2)

figure
I3 = imadjust(I2);
imshow(I3);

figure
bw = imbinarize(I3);
bw = bwareaopen(bw, 50);
imshow(bw)

cc = bwconncomp(bw, 4)