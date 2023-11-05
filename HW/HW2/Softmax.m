function [Pr] = Softmax(a,X_test)
%SOFTMAX Summary of this function goes here
%   Detailed explanation goes here
sigma=@(x)1/(1+exp(-x));
Pr = zeros(400,1);
for i = 1: 400
    Pr(i) = sigma(a(i,:)*X_test');
end
Pr = Pr./sum(Pr);
end

