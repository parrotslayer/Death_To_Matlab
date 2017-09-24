% takes an image and rotates it by either 90, 180 or 270 degrees at random
function output = rot_image(im)

r = randi([1,3],1,1);

output = rot90(im,r);

end