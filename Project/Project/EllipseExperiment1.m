clear
clc
hold on
close all

a = 1 ;
b = 1 ;
A = [1/a^2,...
    0;...
    0,...
    1/b^2];
% Centre is x;y
c =[1;0];
grid on
grid minor
Ellipse_plot(A,c)
x1 = [0;0.5];
x2 = [2;0.5];
x = linspace(x1(1),x2(1),100);
y = linspace(x1(2),x2(2),100);
plot(x,y)
hold on 
collision = CheckCollision(x1,x2,A,c);
if collision == 1
    disp('Collision!')
end
xlim([-5 5])
ylim([-5 5])