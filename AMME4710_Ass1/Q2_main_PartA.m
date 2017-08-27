clc
clear all
close all

%Load validation data
load('legobrick_validation.mat')
load('legobrickjoined_validation.mat')

% Load the image
%lego = imread('legomontage3.png');
lego_num = ['001';'002';'003';'004';'005';'006';'007';'008';'009';'010';...
    '011';'012';'013';'014';'015';'016';'017'];
%choose the image to display
I = 1;

%counters
Yes_Colour = zeros(17,6);
Yes_Center = zeros(17,6);
Yes_Box = zeros(17,6);
True_Pos = zeros(17,6);
True_Pos_Col = zeros(17,6);
False_Neg = zeros(17,1);
False_Pos = zeros(17,1);

%% Begin looping for all images
for I = 1:17

filename = ['legobricks',lego_num(I,:),'.jpg'];
lego = imread(filename);

%Preallocate arrays
% Red, Dark Green, Blue, Light Green, Yellow, Orange
All_centroids = NaN(6,2);
All_boxes = zeros(6,4);

figure
imshow(lego)
title(['Image Loaded: ',filename])

drawbox = 0;

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
[N,trump] = size(regions_Red);
if N > 0
    
    %Store information in arrays
    centroid = regions_Red.Centroid;
    All_centroids(1,:) = centroid;     
    box = regions_Red.BoundingBox;
    All_boxes(1,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Red';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
    
    %count number of boxes drawn
    drawbox = drawbox+1;
end

%% Dark Green Filtering
% Apply filter
[BW,DGreen] = Montage_DGreen2(lego);

%thresholds
min = 400;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Take largest block (no conflicts)
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_DGreen = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a dark green block detected
[N,trump] = size(regions_DGreen);
if N > 0
    
    %Store information in arrays
    centroid = regions_DGreen.Centroid;
    All_centroids(2,:) = centroid;     
    box = regions_DGreen.BoundingBox;
    All_boxes(2,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Dark Green';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
    %count number of boxes drawn
    drawbox = drawbox+1;
end

%% Blue Filtering
% Apply filter
[BW,Blue] = Montage_Blue(lego);

%thresholds
min = 400;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);
% 
% figure
% imshow(Blue)
% regionsos = regionprops(BW, 'Centroid', 'Area', 'BoundingBox')

% Take largest block, high confidence only wanted block is there
BW3 = bwareafilt(BW2,1);

% Details on the final block
regions_Blue = regionprops(BW3, 'Centroid', 'Area', 'BoundingBox');

%check if there is a blue block detected
[N,trump] = size(regions_Blue);
if N > 0
    
    %Store information in arrays
    centroid = regions_Blue.Centroid;
    All_centroids(3,:) = centroid;     
    box = regions_Blue.BoundingBox;
    All_boxes(3,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Blue';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
    %count number of boxes drawn
    drawbox = drawbox+1;
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
[N,trump] = size(regions_LGreen);
if N > 0
    
    %Store information in arrays
    centroid = regions_LGreen.Centroid;
    All_centroids(4,:) = centroid;     
    box = regions_LGreen.BoundingBox;
    All_boxes(4,:) = box; 
    
    %Create bounding box
    rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
    string = 'Light Green';
    text(box(1),box(2)+100,string,'Color','White','FontSize',14)
    %count number of boxes drawn
    drawbox = drawbox+1;
end

%% Yellow Filtering
% Apply filter
[BW,Yellow] = Montage_Yellow1(lego);
%[BW,Yellow] = Montage_Yellow(lego);
%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

% Details on the multiple block
regions_Yellow = regionprops(BW2, 'Centroid', 'Area','BoundingBox');

%check if there is a yellow block detected
[N,trump] = size(regions_Yellow);
if N > 0
    
    for i = 1:N   
        % Find overlapping areas;
        box = regions_Yellow(i).BoundingBox;   
        overlap = 0;
        for k = 1:6
            % find area of overlapping bounding boxes
            area = 0;
            area = rectint(All_boxes(k,:),box);
            %check if any bounding boxes overlap (not with itself)
            if area > 1 && k ~= 5
                %disp('Overlap Found')
                overlap = overlap + 1;
                break
            end
        end        
    end
    
    % dont make bounding box if overlap?
    if overlap < 1 
        %Create bounding box
        rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
        string = 'Yellow';
        text(box(1),box(2)+100,string,'Color','White','FontSize',14)
        % store data in array
        All_boxes(5,:) = box;
        All_centroids(5,:) = regions_Yellow(i).Centroid;
        %count number of boxes drawn
        drawbox = drawbox+1;

    end
    
end

%% Orange Filtering
% Apply filter
[BW,Orange] = MontageOrange(lego);

%thresholds
min = 500;
max = 20000;
BW2 = bwareafilt(BW,[min,max]);

%  figure
%  imshow(BW2)
%  regionsos = regionprops(BW, 'Centroid', 'Area', 'BoundingBox')

% Details on the multiple block
regions_Orange = regionprops(BW2, 'Centroid', 'Area','BoundingBox');

%check if there is a yellow block detected
[N,trump] = size(regions_Orange);
if N > 0
    
    for i = 1:N   
        % Find overlapping areas;
        box = regions_Orange(i).BoundingBox;   
        overlap = 0;
        %check the bounding box against all other boxes
        for k = 1:6
            % find area of overlapping bounding boxes
            area = 0;
            area = rectint(All_boxes(k,:),box);
            %check if any bounding boxes overlap (not with itself)
            if area > 1 && k ~= 6
                disp('Overlap Found')
                overlap = overlap + 1;
                break
            end
        end
        %if passes check is an ok box
        if overlap < 1
            break;
        end
    end
    
    % dont make bounding box if overlap?
    if overlap < 1 
        %Create bounding box
        rectangle('Position',[box(1),(box(2)),box(3),box(4)],'LineWidth',2,'EdgeColor','r');
        string = 'Orange';
        text(box(1),box(2)+100,string,'Color','White','FontSize',14)
        % store data in array
        All_boxes(6,:) = box;
        All_centroids(6,:) = regions_Orange(i).Centroid;
        %count number of boxes drawn
        drawbox = drawbox+1;

    end
    
end

%% validation of data
%get validation data from struct for the specific image
Check_Colour = 0;
Check_Colour = validation_data(I).colours;
Check_Center = validation_data(I).center;
Check_Box = 0;
Check_Box = validation_data(I).box_size;

[trueblocks,Donald] = size(Check_Center);
%% Check if there is a True Block
for i = 1:trueblocks
    %check if we have the right colors
    if strcmp(Check_Colour{i},'red') && All_boxes(1,1) ~= 0
        Yes_Colour(I,1) = 1;
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(1,1) &&...
                Check_Center(i,1) + 50 > All_centroids(1,1) &&...
                Check_Center(i,2) - 50 < All_centroids(1,2) &&...
                Check_Center(i,2) + 50 > All_centroids(1,2)
            
            Yes_Center(I,1) = 1;            
        end        
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(1,3) &&...
                Check_Box(i,1) + 50 > All_boxes(1,3) &&...
                Check_Box(i,2) - 50 < All_boxes(1,4) &&...
                Check_Box(i,2) + 50 > All_boxes(1,4)
            
            Yes_Box(I,1) = 1;
            %All tests passed for true positive
            True_Pos(I,1) = 1;
        end
        
    elseif strcmp(Check_Colour{i},'darkgreen') && All_boxes(2,1) ~= 0
        Yes_Colour(I,2) = 1;        
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(2,1) &&...
                Check_Center(i,1) + 50 > All_centroids(2,1) &&...
                Check_Center(i,2) - 50 < All_centroids(2,2) &&...
                Check_Center(i,2) + 50 > All_centroids(2,2)
            
            Yes_Center(I,2) = 1;             
        end
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(2,3) &&...
                Check_Box(i,1) + 50 > All_boxes(2,3) &&...
                Check_Box(i,2) - 50 < All_boxes(2,4) &&...
                Check_Box(i,2) + 50 > All_boxes(2,4)
            
            Yes_Box(I,2) = 1;
            %All tests passed for true positive
            True_Pos(I,2) = 1;
        end
        
    elseif strcmp(Check_Colour{i},'blue') && All_boxes(3,1) ~= 0
        Yes_Colour(I,3) = 1;
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(3,1) &&...
                Check_Center(i,1) + 50 > All_centroids(3,1) &&...
                Check_Center(i,2) - 50 < All_centroids(3,2) &&...
                Check_Center(i,2) + 50 > All_centroids(3,2)
            
            Yes_Center(I,3) = 1;            
        end
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(3,3) &&...
                Check_Box(i,1) + 50 > All_boxes(3,3) &&...
                Check_Box(i,2) - 50 < All_boxes(3,4) &&...
                Check_Box(i,2) + 50 > All_boxes(3,4)
            
            Yes_Box(I,3) = 1;
            %All tests passed for true positive
            True_Pos(I,3) = 1;
        end
        
    elseif strcmp(Check_Colour{i},'lightgreen') && All_boxes(4,1) ~= 0
        Yes_Colour(I,4) = 1;
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(4,1) &&...
                Check_Center(i,1) + 50 > All_centroids(4,1) &&...
                Check_Center(i,2) - 50 < All_centroids(4,2) &&...
                Check_Center(i,2) + 50 > All_centroids(4,2)
            
            Yes_Center(I,4) = 1;            
        end
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(4,3) &&...
                Check_Box(i,1) + 50 > All_boxes(4,3) &&...
                Check_Box(i,2) - 50 < All_boxes(4,4) &&...
                Check_Box(i,2) + 50 > All_boxes(4,4)
            
            Yes_Box(I,4) = 1;
            %All tests passed for true positive
            True_Pos(I,4) = 1;
        end
        
    elseif strcmp(Check_Colour{i},'yellow') && All_boxes(5,1) ~= 0
        Yes_Colour(I,5) = 1;
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(5,1) &&...
                Check_Center(i,1) + 50 > All_centroids(5,1) &&...
                Check_Center(i,2) - 50 < All_centroids(5,2) &&...
                Check_Center(i,2) + 50 > All_centroids(5,2)
            
            Yes_Center(I,5) = 1;            
        end
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(5,3) &&...
                Check_Box(i,1) + 50 > All_boxes(5,3) &&...
                Check_Box(i,2) - 50 < All_boxes(5,4) &&...
                Check_Box(i,2) + 50 > All_boxes(5,4)
            
            Yes_Box(I,5) = 1;
            %All tests passed for true positive
            True_Pos(I,5) = 1;
        end
        
    elseif strcmp(Check_Colour{i},'orange') && All_boxes(6,1) ~= 0
        Yes_Colour(I,6) = 1;
        %check centroids
        if Check_Center(i,1) - 50 < All_centroids(6,1) &&...
                Check_Center(i,1) + 50 > All_centroids(6,1) &&...
                Check_Center(i,2) - 50 < All_centroids(6,2) &&...
                Check_Center(i,2) + 50 > All_centroids(6,2)
            
            Yes_Center(I,6) = 1;           
        end
        %check bounding box
        if Check_Box(i,1) - 50 < All_boxes(6,3) &&...
                Check_Box(i,1) + 50 > All_boxes(6,3) &&...
                Check_Box(i,2) - 50 < All_boxes(6,4) &&...
                Check_Box(i,2) + 50 > All_boxes(6,4)
            
            Yes_Box(I,6) = 1;
            %All tests passed for true positive
            True_Pos(I,6) = 1;            
        end
    end
end

% Check to see if correct number of bricks detected
if sum(Yes_Colour(I,:)) < trueblocks
    %add number of missed blocks to the counter
    False_Neg(I) = (trueblocks - sum(Yes_Colour(I,:)));
end
%Check to see if too many bricks detected
if drawbox > trueblocks
    %add number of fake blocks to the counter
    False_Pos(I) = drawbox - sum(Yes_Colour(I,:));    
end
%% Check if Wrong Colour
%check for each colour
for k = 1:trueblocks
    if strcmp(Check_Colour{k},'red')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 1
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    
    if strcmp(Check_Colour{k},'darkgreen')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 2
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    
    if strcmp(Check_Colour{k},'blue')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 3
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    
    if strcmp(Check_Colour{k},'lightgreen')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 4
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    if strcmp(Check_Colour{k},'yellow')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 5
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    
    if strcmp(Check_Colour{k},'orange')
        %loop for all blocks
        for j = 1:6
            %check if any center positions match
            if Check_Center(k,1) - 50 < All_centroids(j,1) &&...
                    Check_Center(k,1) + 50 > All_centroids(j,1) &&...
                    Check_Center(k,2) - 50 < All_centroids(j,2) &&...
                    Check_Center(k,2) + 50 > All_centroids(j,2)
                %now check if the index is correct
                if j ~= 6
                    %increment counter
                    True_Pos_Col(I) = True_Pos_Col(I) + 1;
                end
            end
        end
    end
    
end
    
end

%% print out stats

truepos = sum(sum(True_Pos))/17/6*100;
disp(['True Positives = ', num2str(truepos),' %'])

trueposcol = sum(sum(True_Pos_Col))/17/6*100;
disp(['True Positives But Colour Wrong = ', num2str(trueposcol),' %'])

falsepos = sum(sum(False_Pos))/17/6*100;
disp(['False Positives = ', num2str(falsepos),' %'])

falseneg = sum(sum(False_Neg))/17/6*100;
disp(['False Negatives = ', num2str(falseneg),' %'])

leftovers = 100 - truepos - trueposcol - falsepos - falseneg;
disp(['Fully Bricked (leftovers) = ', num2str(leftovers),' %'])