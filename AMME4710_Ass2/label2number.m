% take in a string and converts it to a number for the confusion matrix
function output = label2number(label)

if strcmp(label, 'ball_pit')
    output = 1;
elseif strcmp(label, 'desert')
    output = 2;
elseif strcmp(label, 'park')
    output = 3;
elseif strcmp(label, 'road')
    output = 4;
elseif strcmp(label, 'sky')
    output = 5;
elseif strcmp(label, 'snow')
    output = 6;
elseif strcmp(label, 'urban')
    output = 7;
else
    output = 0;
end

end