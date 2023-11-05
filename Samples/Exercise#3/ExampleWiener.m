clear; clc; close all;

%% General Parameters:
fs = 4e3;   % 4KHz
fc = 120;   % Hz
Tmax = 1;   % seconds (duration of the signal)
s_n = 2;
Ns = Tmax * fs; % signal length


% Define the signals:
t = (0:Ns-1).' /fs;    % time
y = cos(2*pi*fc*t);

% Iteration tests:
Nvec = 5 : 5 : 500;
SNR = -20:20;
Ntest = length(SNR);
Nmc = 100;

MSE = zeros(Ntest, Nmc);

for i = 1 : Ntest
    
    N = 200;
    s_n = 1/db2pow(SNR(i));
    
    for j = 1 : Nmc
        
        x = y + s_n * randn(Ns,1);
        [y_est, h_est, Err_est] = wiener_filter(x, y, N);
        
        MSE(i,j) = mean(abs(Err_est(N+1:end)).^2);
        
    end
    
    
end

hold on;%figure;
plot(SNR, mean(MSE, 2), 'o-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('SNR'); ylabel('MSE');
set(gca, 'FontSize', 24, 'YScale', 'log');

%xlabel('Wiener Order N');

% figure; 
% subplot(3,1,1); hold on;
% plot(t, x, 'k'); plot(t, y, 'r');
% title('Wiener Filter Problem');
% legend({'Noisy Signal', 'Desired'});
% 
% subplot(3,1,2); hold on;
% plot(t, y_est, 'k'); plot(t, y, 'r');
% title('Wiener Filter Problem');
% legend({'Estimated', 'Desired'});
% 
% subplot(3,1,3); hold on;
% plot(t, Err_est, 'k');
% title('Wiener Filter Problem');
% legend({'Residual'});



%% Usefull functions:

function [y_estimated, h, e_estimated] = wiener_filter(x, y, N)

    Xw = fft(x(1:N)) / N;
    Yw = fft(y(1:N)) / N;
    
    Rxx = N * real(ifft(Xw .* conj(Xw)));
    Rxy = N * real(ifft(Xw .* conj(Yw)));
    
    Rxx = toeplitz(Rxx);
    Rxy = Rxy';
    
    h = Rxy / Rxx;

    y_estimated = fftfilt(h,x);
    
    e_estimated = y - y_estimated;
end