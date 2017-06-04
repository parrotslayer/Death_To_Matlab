% KNN example 2
clear all
clc

point_coordinates =   [11.4179  103.1400
   16.7710   10.6691;
   16.6068  119.7024;
   25.1379   74.3382;
   30.3651   23.2635;
   31.7231  105.9109;
   31.8653   36.9388];

[idx,d] = knnsearch(point_coordinates, point_coordinates, 'k', 2);
idx = idx(:,2)
d = d(:,2)