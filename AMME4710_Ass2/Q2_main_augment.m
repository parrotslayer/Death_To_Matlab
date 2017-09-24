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
tbl = countEachLabel(imds);
%%
% Because |imds| above contains an unequal number of images per category,
% let's first adjust it, so that the number of images in the training set
% is balanced.

minSetCount = min(tbl{:,2}); % determine the smallest amount of images in a category;

% Use splitEachLabel method to trim the set.
imds_train = splitEachLabel(imds, minSetCount, 'randomize');

% Notice that each set now has exactly the same number of images.
countEachLabel(imds);

% Divide data into training and testing set
fold_ratio = 0.3;

fold = 5;
%for fold = 5:5
fold_ratio = fold/10;
[trainingSet, testSet] = splitEachLabel(imds, fold_ratio, 'randomize');

%% Make the machine learning model or just load it
% Make the model

make_model = 1;
inc = 1;
featureVector = NaN(256*fold_ratio*2,103917);
featureVector2 = NaN(256*(1-fold_ratio),103917);
if make_model == 1

    for aug = 1:2
        for i = 1:length(trainingSet.Files)
            %for i = 1:5
            %extract features using extractHOG
            clear image
            clear featureVector_col
            clear featureVector_Hist
            % Augment image by random rotation or crop
            if aug == 1
                image = trainingSet.readimage(i);
            elseif aug == 2
                image1 = trainingSet.readimage(i);
                image = rot_image(image1);
            elseif aug == 3
                image1 = trainingSet.readimage(i);
                image = crop_image(image1);
            end

            [featureVector_HOG,hogVisualization] = extractHOGFeatures(image);
            %normalise the feature vector
            featureVector_HOG = featureVector_HOG./max(featureVector_HOG);

            [rows,cols,pages] = size(image);
        % Get raw pixel data of colour
        k = 1;
        KS = 0;
        for r = 1:rows
            for c = 1:cols
                for p = 1:pages
                    %subsample the image data by a factor KS
                    if KS == 4
                        %make 1D array of the colour image
                        featureVector_col(k) = image(r,c,p);
                        k = k + 1;
                        KS = 0;
                    else
                        KS = KS + 1;
                    end                    
                end
            end
        end
            %normalise the featureVector and convert to double
            featureVector_col = double(featureVector_col./max(featureVector_col));

            % Use historgram data over the RGB channels. Could do over HSV if
            % wanted to.
            % Repeat for R, G, B channels with X bins each
            num_bins = 10000;
            [counts_R,bins] = imhist(image(:,:,1),num_bins);
            [counts_G,bins] = imhist(image(:,:,2),num_bins);
            [counts_B,bins] = imhist(image(:,:,3),num_bins);
            featureVector_Hist = [counts_R',counts_G',counts_B'];

            %normalise the featureVector along ALL the 3 channels.
            % Could normalise for each channel seperately?
            featureVector_Hist = double(featureVector_Hist./max(featureVector_Hist));

            % Concatenate the different feature vectors
            featureVector(inc,:) = [featureVector_HOG, featureVector_col,featureVector_Hist];
            inc = inc + 1;
        end

    end
    training_labels_aug = [trainingSet.Labels;trainingSet.Labels];
    % Train an ECOC multiclass model using the default options. SVM
    Mdl = fitcecoc(featureVector,training_labels_aug);
    % save the generated model
    save('Mdl_HOG_col_hist.mat','Mdl');
else
    load('Mdl_HOG_col.mat')
end

%% Pass in the testing data
% take first remaining folds and classify the data
for i = 1:length(testSet.Files)
    %extract features using extractHOG
    clear image
    clear featureVector_col
    clear featureVector_Hist
    image = testSet.readimage(i);
    [featureVector_HOG,hogVisualization] = extractHOGFeatures(image);
    %normalise the feature vector
    featureVector_HOG = featureVector_HOG./max(featureVector_HOG);

    [rows,cols,pages] = size(image);
    % Get raw pixel data of colour
    k = 1;
    KS = 0;
    for r = 1:rows
        for c = 1:cols
            for p = 1:pages
                %subsample the image data by a factor KS
                if KS == 4
                    %make 1D array of the colour image
                    featureVector_col(k) = image(r,c,p);
                    k = k + 1;
                    KS = 0;
                else
                    KS = KS + 1;
                end                    
            end
        end
    end

    %normalise the featureVector and convert to double
    featureVector_col = double(featureVector_col./max(featureVector_col));

    % Use historgram data over the RGB channels. Could do over HSV if
    % wanted to.
    % Repeat for R, G, B channels with X bins each
    num_bins = 10000;
    [counts_R,bins] = imhist(image(:,:,1),num_bins);
    [counts_G,bins] = imhist(image(:,:,2),num_bins);
    [counts_B,bins] = imhist(image(:,:,3),num_bins);
    featureVector_Hist = [counts_R',counts_G',counts_B'];

    %normalise the featureVector along ALL the 3 channels.
    % Could normalise for each channel seperately?
    featureVector_Hist = double(featureVector_Hist./max(featureVector_Hist));

    % Concatenate the different feature vectors
    featureVector2(i,:) = [featureVector_HOG, featureVector_col,featureVector_Hist];

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
        i_pred = label2number(char(label_pred(i)));
        i_true = label2number(char(testSet.Labels(i)));
        conf(i_pred,i_true) = conf(i_pred,i_true) + 1;
    else
        %figure
        %imshow(testSet.readimage(i))
        %title(['P = ',char(label_pred(i)),' T = ',char(testSet.Labels(i))])
        %increment confusion matrix
        i_pred = label2number(char(label_pred(i)));
        i_true = label2number(char(testSet.Labels(i)));
        conf(i_pred,i_true) = conf(i_pred,i_true) + 1;
    end
end

%% Show evaluation
disp(['Fold Ratio = ',num2str(fold_ratio*100),' %']);
true_pos_rate = true_pos/length(testSet.Files)*100;
disp(['True Positive Rate = ', num2str(true_pos_rate),' %'])

for i = 1:7
    %calc precision (across)
    precision(i) = conf(i,i)/sum(conf(i,:));
    %calc recall (down)
    recall(i) = conf(i,i)/sum(conf(:,1));
end
F1 = 2*(precision.*recall)/(precision+recall)*100;
disp(['F1 Score = ', num2str(F1),' %'])
disp('Confusion Matrix')
disp(conf)
results = zeros(8);
results(2:8,1) = precision;
results(1,2:8) = recall;
results(2:8,2:8) = conf;
    
%end    %end fold for loop
%% show images
% i = 10;
% for i = 100:110
%     figure
%     imshow(testSet.readimage(i))
%     title(['P = ',char(label_pred(i)),' T = ',char(testSet.Labels(i))])
% end