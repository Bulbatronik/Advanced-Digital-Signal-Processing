function [H] = Class(point)
%CLASS Summary of this function goes here
%   Return the class
delta = 0.5;
x = point(1)+5;
y = point(2) - 20;


x1 = floor(x/delta);
y1 = ceil(y/delta);

a = x1;
b = 21-y1;

H = 20* a +b;
end

