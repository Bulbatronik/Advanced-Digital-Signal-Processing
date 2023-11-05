%% WS#2 - Solutions
clc, clear, close all;

% E1.1: Noptimal = 1, for both Wiener and LMS, MU = 0.015 for LMS
% E1.2: Noptimal = 6 for Wiener, Noptimal = 7, MU = 0.012 for LMS
% E1.3: Noptimal = 10 for Wiener, Noptimal = 5, MU = 0.006 for LMS

%% E1:
load('Sin_3.mat');

Nsig = length(x3);
Nsig_t = round(0.2 * Nsig);
x_t = x3(1:Nsig_t);
n_t = n3(1:Nsig_t);

Nvect = 1 : 25;
Ntest = length(Nvect);

mu_vect = 0.001:0.001:0.015;
Nmu = length(mu_vect);

MSE_vs_N_wiener = zeros(Ntest, 1);
MSE_vs_N_lms = zeros(Ntest, Nmu);

for n = 1 : Ntest
   
    N = Nvect(n);
    
    % Wiener solution:
    [y_estimated_w, h_estimated_w, e_estimated_w] = wiener_filter(x_t, n_t, N);
    MSE_vs_N_wiener(n) = mean(e_estimated_w(N+1:end).^2);
    
    % LMS solution:
    for m = 1 : Nmu
        mu = mu_vect(m);
        ho = zeros(N,1);
        [y_estimated_lms, h_estimated_lms, e_estimated_lms] = LeastMeanSquare(x_t, n_t, mu, N, ho);
        MSE_vs_N_lms(n, m) = nanmean(e_estimated_lms(N+1:end).^2);
    end
    
end

%% Plots:
figure, hold on, grid on, box on;
plot(Nvect, MSE_vs_N_wiener, '-o', 'LineWidth', 2)
xlabel('Filter Length (N)'), ylabel('MSE')
set(gca, 'FontSize', 20)

figure; hold on;
plot(Nvect, MSE_vs_N_lms, '-x', 'LineWidth', 2)
xlabel('Filter Length (N)'), ylabel('MSE')
set(gca, 'FontSize', 20)

%% E2
clc; clear; close all;

load('Sdoppler.mat');

Nsig = length(sdoppler);
Nseg = 500;
Nwin = round(Nsig / Nseg);


Vsound = 335;               % [m/s] acousting wave speed over the air

speedRange = (30:50) / 3.6; % [m/s] speed range in urban scenario


Nrange = length(speedRange);

Pinitial = [100, 5];
Pfinal = [-100, 5];
dist = pdist2(Pinitial, Pfinal);
ds = [dist / Nwin, 0];

iS = 1;

Ppast = Pinitial;
Pdest = [0 0];

EstSpeed = [];

for ns = 1 : Nwin-1
    
    Pcur = Ppast - ds;
    Acur = atan2(Pcur(2) - Pdest(2), Pcur(1) - Pdest(1));
    
    str_seg = strue((1:Nseg) + (ns-1) * Nseg);
    
    RMax = zeros(Nrange, 1);
    
    iF = zeros(Nrange, 1);
    
    for nr = 1 : Nrange
        fdop = Vsound / (Vsound + (speedRange(nr) * cos(Acur)));
        
        t_true = 1 : Nseg;
        t_interp = linspace(1, Nseg, Nseg * fdop);
        
        str_seg_dop = interp1(t_true, str_seg, t_interp);
        
        % The length of the interpolated signal change depending on the
        % doppler shift, it is important to span the sdoppler and the strue
        % taking into account the different lengths
        iF(nr) = length(str_seg_dop);
        
        sdop_seg = sdoppler(iS:iS+iF(nr)-1);
        
        r_st_sd = xcorr(str_seg_dop, sdop_seg);
        
        RMax(nr) = r_st_sd(iF(nr));

    end
    
    % This was the error, by not considering the proper span og the signal,
    % you may correlate two different parts of the signal which experience
    % different dopplers. Below, I span the sdoppler based on the length
    % that is defined by the previously estimated length:
    iS = iS + round(mean(iF(RMax == max(RMax)))); 
    
    % Since I may have multiple speeds with the same corss-correlation, I collect
    % all the speeds which maximeze the corss-correlation.
    EstSpeed = [EstSpeed 3.6*speedRange(RMax == max(RMax))];
    
    Ppast = Pcur;
end

% plots:
histogram(EstSpeed, 20, 'Normalization', 'pdf');     % Here, the histogram of estimated speeds
xlabel('Speed [Km/h]'); ylabel('PDF')
set(gca, 'FontSize', 20);

fprintf('The estimated speed is: %.2f [Km/h]\n', mean(EstSpeed));

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