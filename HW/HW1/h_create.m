function [h] = h_create(m,n,K,alpha,b)
%H Summary of this function goes here
%   Creatin the matrix h for a given alpha
h = zeros(m, n, K);
for i= 1:m
    for l = 1:n
        for k = 1:K
            h(i,l,k) = alpha^abs(i-l)*b^k;
        end
    end
end

end

