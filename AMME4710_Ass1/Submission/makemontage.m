
clear all
close all
clc

lego(:,:,:,1) = imread('legobricks001.jpg');
lego(:,:,:,2) = imread('legobricks002.jpg');
lego(:,:,:,3) = imread('legobricks003.jpg');
lego(:,:,:,4) = imread('legobricks004.jpg');
lego(:,:,:,5) = imread('legobricks005.jpg');
lego(:,:,:,6) = imread('legobricks006.jpg');
lego(:,:,:,7) = imread('legobricks007.jpg');
lego(:,:,:,8) = imread('legobricks008.jpg');
lego(:,:,:,9) = imread('legobricks009.jpg');
lego(:,:,:,10) = imread('legobricks010.jpg');
lego(:,:,:,11) = imread('legobricks012.jpg');
lego(:,:,:,12) = imread('legobricks013.jpg');
lego(:,:,:,13) = imread('legobricks014.jpg');
lego(:,:,:,14) = imread('legobricks015.jpg');
lego(:,:,:,15) = imread('legobricks016.jpg');
lego(:,:,:,16) = imread('legobricks017.jpg');
%lego(:,:,:,17) = imread('legobricks017.jpg');

% lego(:,:,:,11) = imread('legobricks011.jpg');
% lego(:,:,:,12) = imread('legobricks012.jpg');
% lego(:,:,:,13) = imread('legobricks013.jpg');
% lego(:,:,:,14) = imread('legobricks014.jpg');
% lego(:,:,:,15) = imread('legobricks015.jpg');
% lego(:,:,:,16) = imread('legobricks016.jpg');
% lego(:,:,:,17) = imread('legobricks017.jpg');

%% 
montage(lego)
