function hit = CheckCollision(vec1,vec2,A,c)
%Students to add code that returns:
% hit = 1 if the line from x1 to x2 hits the ellipse (x-c)'A(x-c)<=1.
% hit = 0 otherwise (i.e. path is safe)

res = 100; %number of points to divide up the line
hit = 0;

x1 = vec1(1);
y1 = vec1(2);
x2 = vec2(1);
y2 = vec2(2);

%generate res number of evenly spaced points between the two vectors
X = linspace(x1,x2,res);
Y = linspace(y1,y2,res);

for k = 1:res
    X_vec = [X(k); Y(k)];
    test = (X_vec-c)'*A*(X_vec-c);
    
    %test if point is within the ellipse
    if test <= 1
        %if within ellipse, break
        hit = 1;
        break
    end    
end
end


