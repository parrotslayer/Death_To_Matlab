% takes an images and crops into one of 4 quadrants
function output = crop_image(im)

r = randi([1,4],1,1);

[rows, cols, pages] = size(im);

if r == 1
    output = im(1:rows/2, 1:cols/2, :);
elseif r == 2
    output = im(1:rows/2, cols/2:cols, :);
elseif r == 3
    output = im(rows/2:rows, 1:cols/2, :);
elseif r == 4
    output = im(rows/2:rows, cols/2:cols, :);
end

output = imresize(output, [rows cols]);

end