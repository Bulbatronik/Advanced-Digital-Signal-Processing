function [C] = Cov_True(rho, m)
%COV_TRUE Summary of this function goes here
%  For given rho and given m compute rhe true cov matrix
for i = 1:m
    for j = 1 : m
       if i == j
          C(i,j) = 1;
       else 
          C(i,j) = rho;
       end
    end
end

