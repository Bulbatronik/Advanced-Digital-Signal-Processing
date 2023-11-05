%% ADSP 2021 - 2022 / WS #1 - solution

clc; clear; close all;

% Load data:
load('data_1.mat');
load('true_signal.mat');

%% 1.a) given the model: x[n] = f(s[n]) + w[n], estimate the distortion f(.)

P_vec = 1 : 20;             % I'm defining the set of possible polynomial orders

A_est = zeros(P_vec(end));
mse = zeros(P_vec(end), 1);
S = [];

for p = 1 : length(P_vec)
   
    S = [S, s.^p];      % for every loop, I'm considering a higher order of polynomial distortion
    a_est = S \ x1;
    A_est(p,1:length(a_est)) = a_est;
    mse(p) = mean(abs(x1 - S * a_est).^2);
    
end

figure; grid on;
plot(P_vec, mse, '-x', 'LineWidth', 2);
xlabel('Polynomial Order [p]')
ylabel('MSE')
set(gca, 'FontSize', 24, 'YScale', 'log')
ylim([1e-3 1])

%% 1.b) given the model: x[n] = h[n] * s[n] + w[n], estimate the system h[n]

% The size of the entire data set is high, we need to use a truncated
% version. However, the truncation introduce an error and comes with a
% tradeoff between good and fast estimate. To compensate for the error
% due to truncation, we need to use the overlapp and add method [see 
% literature for references]

clc; clear; close all;

% Load data:
load('data_2.mat');
load('true_signal.mat');

Ns = size(s, 1);                % size of the signal
L_vec = 100:100:2000;           % assuming the h[n] response has a length between 100 and 2000
h_est_v = zeros(L_vec(end));

for i = 1 : length(L_vec)
    
    N = L_vec(i);                       % considered filter length
    T_oa = 8 * N;                       % I'm truncating the signal x to be at least 8x the length of the filter.
    Nstep = T_oa - N + 1;               % length of the truncation
    Nsegments = round(Ns / Nstep) - 1;  % number of possible truncations consdering the entire signal
    H_est = zeros(N,Nsegments-3);

    for n = 2 : Nsegments-2

        i1 = Nstep * n;
        i2 = 1:Nstep;
        i3 = 1:T_oa;

        S = zeros(T_oa, N);
        for ig = 1 : N
            S(ig:ig+Nstep-1, ig) = s(i1+i2);
        end

        s_p1 = s(i1+i2+Nstep);
        S_p1 = zeros(T_oa, N);
        for ig = 1 : N
            S_p1(Nstep+ig:end, ig) = s_p1(1:N-ig);
        end

        h_tmp = (S+S_p1) \ x2(i1+i3);

        H_est(:,n-1) = h_tmp; 
    end
    
    h_est = mean(H_est, 2);
    h_est_v(i,1:length(h_est)) = h_est;
    
    x2_predicted = conv(s, h_est, 'same');
    mse(i) = mean(abs(x2 - x2_predicted).^2);
    
        
end

figure; grid on;
plot(L_vec, mse, '-x', 'LineWidth', 2);
xlabel('System length [L]')
ylabel('MSE')
set(gca, 'FontSize', 24, 'YScale', 'log')

figure;
stem(h_est_v(10,1:1000))

