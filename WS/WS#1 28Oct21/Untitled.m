%% 
clc; clear all;
load data_1.mat
load true_signal.mat
%sound(x1, Fs)
sigw_2 = -20;% dB
w=sigw_2*randn(size(x1));
%% 
p = 5;
S = zeros(length(x1), p);
a = zeros(p, p);
MSE = zeros(p, 1);
%% COMPUTING S
for i =1:length(x1)
    for j = 1:p
        S(i,j) = s(i)^j;
    end
end
%% a = INV(S'*S)*S'*X1
for i= 1:p
    for j = 1:i
        %theta(i,j) =  inv(S(i,j)'* S(i,j))* S(i,j)'*x1;
        mat = pinv(S(:,1:i))*x1;
        a(i,j) = mat(j);
    end
end
%% ESTIMATING SIGNAL
%ELEMENT IN POWER 1 CORRESPONDS TO THE FIRST ELEMENT IN EVERY ROW
s_rec = zeros(length(s), p);
for i = 1:p
    for j = 1:length(s)
        %mat = pinv(a(i,1:i))*x1';
        mat = pinv(a(i,1))*x1';
        s_rec(i,j) = mat(j);
    end
end
%% MSE
for i = 1:p
    MSE(i) = mean((s_rec(:,i)-s).^2);
end
plot(1:p, MSE)