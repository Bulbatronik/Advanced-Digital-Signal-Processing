clear; clc; close all;

%% Parameters:
K = 500;
s_n = 0.04;
t = 0:K-1;


% true filter:
N_true = 4;
h_true = [1 .5 .1 .01]';

% signal generation:
x = randn(K,1);
s_x = var(x);

% noise generation:
n = s_n * randn(K,1);

% Measurements:
y = zeros(K,1);
X = zeros(N_true,1);
for k = 1: K
    X = [x(k,1); X(1:(N_true-1),1)];
    y(k) = h_true' * X(:,1) + n(k);
end


% Test:
N = 100;
mu = 0.04;
ho = randn(N,1);
[x_est, h_est, Err_est] = LeastMeanSquare(y, x, mu, N, ho);


figure; 
subplot(3,1,1); hold on;
plot(t, y, 'k'); plot(t, x, 'r');
title('LMS Filter Problem');
legend({'Noisy Signal', 'Desired'});

subplot(3,1,2); hold on;
plot(t, x_est, 'k'); plot(t, x, 'r');
title('LMS Filter Problem'); ylim([-10 10])
legend({'Estimated', 'Desired'});

subplot(3,1,3); hold on;
plot(t, Err_est, 'k'); 
title('LMS Filter Problem'); ylim([-10 10])
legend({'Residual'});


%% Usefull functions:

function [x_estimated, h_estimated, e_estimated] = LeastMeanSquare(y, x, mu, N, ho)

    K = length(x);
    e_estimated = zeros(K, 1);
    x_estimated = zeros(K, 1);
    h_estimated = zeros(N, K);
    
    h_estimated(:,1) = ho;
    
    prefixedInput = [zeros(N-1,1); y];
    
    for it = 1 : K
       
        regressor = prefixedInput(it+N-1:-1:it,1);
        
        x_estimated(it, 1) = h_estimated(:,it)' * regressor;
        
        e_estimated(it, 1) = x(it) - x_estimated(it, 1);
        
        h_estimated(:,it+1) = h_estimated(:,it) + (mu * conj(e_estimated(it, 1)) * regressor);
        
    end
    
end