% take in a string and converts it to a number for the confusion matrix
function output = label2number(label)

if label == 'ball_pit'
    output = 1;
elseif label == 'desert'
    output = 2;
elseif label == 'park'
    output = 3;
elseif label == 'road'
    output = 4;
elseif label == 'sky'
    output = 5;
elseif label == 'snow'
    output = 6;
elseif label == 'urban'
    output = 7;
else
    output = 0;
end

end