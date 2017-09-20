clc
close all
clear all

%% Load Images
% Instead of operating on all of Caltech 101, which is time consuming, use
% three of the categories: airplanes, ferry, and laptop. The image category
% classifier will be trained to distinguish amongst these six categories.

%rootFolder = fullfile(outputFolder, 'dataset');
%categories = {'airplanes', 'ferry', 'laptop'};

rootFolder = 'C:\Users\Bill\Documents\GITHUB\Death_To_Matlab\AMME4710_Ass2\dataset';

%%
% Create an |ImageDatastore| to help you manage the data. Because
% |ImageDatastore| operates on image file locations, images are not loaded
% into memory until read, making it efficient for use with large image
% collections.
imds = imageDatastore(rootFolder,'IncludeSubfolders',true, 'LabelSource', 'foldernames');

%%
% The |imds| variable now contains the images and the category labels
% associated with each image. The labels are automatically assigned from
% the folder names of the image files. Use |countEachLabel| to summarize
% the number of images per category.
tbl = countEachLabel(imds)
%%
% Because |imds| above contains an unequal number of images per category,
% let's first adjust it, so that the number of images in the training set
% is balanced.

minSetCount = min(tbl{:,2}); % determine the smallest amount of images in a category;

% Use splitEachLabel method to trim the set.
imds_train = splitEachLabel(imds, minSetCount, 'randomize');

% Notice that each set now has exactly the same number of images.
countEachLabel(imds)

% Divide data into training and testing set
[trainingSet, testSet] = splitEachLabel(imds, 0.3, 'randomize');

%% Tute 7
% Make the model
%featureVector = NaN(length(trainingSet.Files),250000);

for i = 1:length(trainingSet.Files)
    %for i = 1:5
    %extract features using extractHOG
    clear image
    image = trainingSet.readimage(i);
    [featureVector_HOG,hogVisualization] = extractHOGFeatures(image);
    %normalise the feature vector
    featureVector_HOG = featureVector_HOG./max(featureVector_HOG);
    
    [rows,cols,pages] = size(image);
    % Get raw pixel data of colour
    k = 1;
    for r = 1:rows
        for c = 1:cols
            for p = 1:pages
                %make 1D array of the colour image
                featureVector_col(k) = image(r,c,p);
                k = k + 1;
            end
        end
    end
    
    %normalise the featureVector and convert to double
    featureVector_col = double(featureVector_col./max(featureVector_col));
    
    % Concatenate the different feature vectors
    featureVector(i,:) = [featureVector_HOG, featureVector_col];
    
    %convert featureVector to type double because gets angry
    %disp(i)
    %     figure
    %     subplot(1,2,1);
    %     imshow(image);
    %     subplot(1,2,2);
    %     plot(hogVisualization);
end


% Train an ECOC multiclass model using the default options. SVM
Mdl = fitcecoc(featureVector,trainingSet.Labels);
save('Mdl_HOG_col.mat','Mdl');

%% Pass in the testing data
% take first remaining folds and classify the data
for i = 1:length(testSet.Files)
    %extract features using extractHOG
    clear image
    image = testSet.readimage(i);
    [featureVector_HOG,hogVisualization] = extractHOGFeatures(image);
    %normalise the feature vector
    featureVector_HOG = featureVector_HOG./max(featureVector_HOG);
    
    [rows,cols,pages] = size(image);
    % Get raw pixel data of colour
    k = 1;
    for r = 1:rows
        for c = 1:cols
            for p = 1:pages
                %make 1D array of the colour image
                featureVector_col(k) = image(r,c,p);
                k = k + 1;
            end
        end
    end
    
    %normalise the featureVector and convert to double
    featureVector_col = double(featureVector_col./max(featureVector_col));
    
    % Concatenate the different feature vectors
    featureVector2(i,:) = [featureVector_HOG, featureVector_col];
    
end

% Pass features into predict. Returns vector with predicted
label_pred = predict(Mdl,featureVector2);

%Assess whether true positive or not
true_pos = 0;
conf = zeros(7,7);
for i = 1:length(testSet.Files)
    if label_pred(i) == testSet.Labels(i)
        true_pos = true_pos + 1;
        %increment confusion matrix
        %conf(label_pred(i)+1,testSet.Labels(i)+1) = conf(label_pred(i)+1,testSet.Labels(i)+1) + 1;
        %plot incorrect ones
    else
        %         figure
        %         imshow(im_sub2(:,:,i))
        %         title(['P = ',num2str(label_pred(i)),'  T = ',num2str(labels_sub2(i))])
        %increment confusion matrix
        %conf(label_pred(i)+1,testSet.Labels(i)+1) = conf(label_pred(i)+1,testSet.Labels(i)+1) + 1;
    end
end

%% Show evaluation
true_pos_rate = true_pos/length(testSet.Files)*100;
disp(['True Positive Rate = ', num2str(true_pos_rate),' %'])
disp('Psyduck used Confusion Matrix')
disp(conf)

%% show images
i = 10;
for i = 100:110
    figure
    imshow(testSet.readimage(i))
    %title([label_pred(i)])
    disp(label_pred(i))
    %disp(testSet.Labels(i))
end