%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     I = ro*N*s    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% where:
% ro -> albedo (reflectance) - 
% N -> Normal
% s -> direction of the light

% im_array - has all the images

%only first pixel for now


% I = s*g
% b = Ax    ...  x = A^-1 b  or  A^-1\b

load('facedata_yaleB07.mat')

[h,w,N] = size(im_array);

A = light_dirs;

for i = 1:h
    
    for j = 1:w
        
     
    b = reshape(im_array(i,j,:),[N,1]);  % 1x1x64 -> 64x1
    x = A\(double(b )/255);
    
    rho(i,j) = norm(x);
    nhat = x/rho(i,j);
    New(i,j,1) = nhat(1);
    New(i,j,2) = nhat(2);
    New(i,j,3) = nhat(3);
    
    end
    
end

figure
imagesc(rho)
colormap(gray)