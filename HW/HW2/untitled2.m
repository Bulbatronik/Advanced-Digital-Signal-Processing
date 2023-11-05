%point = [-4.9, 29.9];%1
%point = [+4.9, 29.9];%381
%point = [+4.9, 20.1];%400
delta = 0.5;
x = point(1)+5;
y = point(2) - 20;


x1 = floor(x/delta);
y1 = ceil(y/delta);

a = x1;
b = 21-y1;

area = 20* a +b;
% if point(1) > 0
%     x =  floor(2*(point(1)+5));%/delta);
% else
%     x =  ceil(2*(point(1)+5));%/delta);
% end
% 
% if point(2) > 0
%     y = floor(2*(point(2) -20));%/delta)-20;
% else
%     y =  ceil(2*(point(2)- 20));%/delta)-20;
% end

% x = point(1)+5;
% y = point(2) - 20;
% 
% x  = round(x * 2);
% y = round(y * 2);

% a = abs(x);
% b = abs(y);
% 
% area = 20* a + b ;
% if b ==0
%     b=1;
% end

% 
% 
% a = abs(-10 - x);
% b = abs(60 - y);
% 
% area = (a-1)*20 + b;
