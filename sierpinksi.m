% Sierpinski fractal attempt.
clc
clear
% Initialise an equilateral triangle 
vertex1 = [0 0];
vertex2 = [5 0];
vertex3 = [2.5 4.33];
triangle = [vertex1;vertex2;vertex3];
scatter(triangle(:,1),triangle(:,2),1)
% Initialise first point and fractal db
point = [2.5 0];
fractal = point;
tic
for i = 1:1000000
    rand = randsample(1:3,1);
    if rand == 1
        vertex = vertex1;
    elseif rand == 2
        vertex = vertex2;
    elseif rand == 3
        vertex = vertex3;
    end
    mid_dir = [(point(1)-vertex(1))/2 (point(2)-vertex(2))/2];
    new_point = point - mid_dir;
    fractal = [fractal;new_point];
    point = new_point;
end
toc
hold on
scatter(fractal(:,1),fractal(:,2),1)