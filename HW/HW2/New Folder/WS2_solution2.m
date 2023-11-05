clc; clear; close all;

load('data_2.mat')

N = 10;
mu = 0.001;
lam = 0.99;
[y_lms, h_lms, e_lms] = LMS(w1, x1, mu, N, zeros(N,1));
[y_rls, h_rls, e_rls] = EW_RLS(w1, x1, lam, N, zeros(N,1));

figure; hold on;
plot(abs(e_lms(:,1)).^2, 'LineWidth', 2)
plot(abs(e_rls(:,1)).^2, 'LineWidth', 2)
xlabel('Discrete-Time index [n]')
ylabel('\xi^2[n]')
set(gca, 'FontSize', 24)

%% Useful functions

function [theta_est, h_est, e_est] = LMS(x, theta, mu, N, ho)

    [K, M] = size(theta);
    e_est = zeros(K, M);
    theta_est = zeros(K, M);
    h_est = zeros(N, M, K);
    for m = 1 : M
        h_est(:,m,1) = ho;
    end
    
    prefixedInput = [zeros(N-1,M); x];
    
    for it = 1 : K
       
        for m = 1 : M
            regressor = prefixedInput(it+N-1:-1:it,m);

            ho = h_est(:,m,it);
            theta_est(it, m) = ho' * regressor;

            e_est(it, m) = theta(it,m) - theta_est(it, m);

            ho = ho + (mu * conj(e_est(it, m)) * regressor);
            h_est(:,m,it+1) = ho;
        end
    end
    
end

function [x_est, h_est, e_est] = EW_RLS(x, theta, lambda, N, ho)

    [K, M] = size(x);
    e_est = zeros(K, M);
    x_est = zeros(K, M);
    h_est = zeros(N, M, K);
    for m = 1 : M
        h_est(:,m,1) = ho;
    end
    
    P = eye(N);
    
    prefixedInput = [zeros(N-1,M); x];
    
    for it = 1 : K
       
        for m = 1 : M
            regressor = prefixedInput(it+N-1:-1:it,m);

            ho = h_est(:,m,it);
            x_est(it, m) = ho' * regressor;

            e_est(it, m) = theta(it,m) - x_est(it, m);
            
            g = P * regressor / (lambda + regressor' * P * regressor);
            P = P/lambda - g * regressor' * (1/lambda) * P;
            ho = ho + e_est(it, m) * g;
            h_est(:,m,it+1) = ho;
        end
    end
    
end