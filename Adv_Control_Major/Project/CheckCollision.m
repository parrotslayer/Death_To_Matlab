function hit = CheckCollision(x1,x2,A,c)
gridpoints = 100;
x = linspace(x1(1),x2(1),gridpoints);
y = linspace(x1(2),x2(2),gridpoints);
linepoints = [x;y];
value = NaN(1,gridpoints);
hit = 0;
for n = 1:length(x)
    value(n) = (linepoints(:,n) - c)' * A * (linepoints(:,n) - c);
    if value(n) <= 1
        hit = 1;
        break
    end
end
end

    
% hit = 1 if the line from x1 to x2 hits the ellipse (x-c)'A(x-c)<=1.
% hit = 0 otherwise (i.e. path is safe)


