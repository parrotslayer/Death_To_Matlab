function hit = CheckCollisionPoint(vec1,A,c)
%Students to add code that returns:
% hit = 1 if the line from x1 to x2 hits the ellipse (x-c)'A(x-c)<=1.
% hit = 0 otherwise (i.e. path is safe)

x1 = vec1(1);
y1 = vec1(2);

X_vec = [x1; y1];
test = (X_vec-c)'*A*(X_vec-c);

%test if point is within the ellipse
if test <= 1
    %if within ellipse, break
    hit = 1;
else
    hit = 0;
end



