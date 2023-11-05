%% 1.3 MIMO deconvolution
clear all; clc;
M = 4;
N = 4;

SNR = -10:2:30; %dB
P = 200;
sigma_X2 = 10.^(SNR/10);

Q = floor(linspace(5, P-1, 20));%Q_min>=K(M-1)+1

rho = 0.1;
alpha = 0.5;

h = zeros(M, N);
for l= 1:N
    for i = 1:M
        h(i,l) = alpha^abs(i-l);
    end
end

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/199'])
    %X = zeros(floor(Q(Q_index)),M);
    for snr = 1:length(sigma_X2)
        for p =1 :50%MC
            x = sigma_X2(snr).* randn(floor(Q(Q_index)),N);
            C = Cov_True(rho, M);%4X4
            w = chol(C, 'lower')* randn(M, floor(Q(Q_index)));%%%%%%%%%%%%%%%%%%%%
           

            y = zeros(floor(Q(Q_index)+1-1), M);
            for j=1:M
                y(:,j) = x*h(j,:)' + w(j,:)';
                
                %X_ml = (X.' / C *X) \ X.' / C;
                crb_h = (x'*x)^-1;
                h_est(j,:) = (crb_h*x'*y(:,j))';%(X_ml*y(:,j))';

                MSE_h(Q_index , snr, p) = mean((h_est(j,:) - h(j,:)).^2, "all");%,"all"
            end
            
             x2 = sigma_X2(snr).* randn(floor(P-Q(Q_index)), N);
             w2 = chol(C, 'lower')* randn(M, floor(P-Q(Q_index)));
              
             y2 = zeros(floor(P-Q(Q_index)), M);
             for i = 1:floor(P-Q(Q_index))
                 y2(i,:) = (h*x2(i,:)' + w2(:,i))';

                 %%H_ml = (h_est.'*h_est) \ h_est.';%uncorrelated between each other
                 
                 %H_ml = (h_est.' /Cov_True(rho, M)*h_est) \ h_est.'* Cov_True(rho, M);
                 crb_x = (h_est' * C^-1 * h_est)^-1;%* h_est.'* Cov_True(rho, M);
                 X_est_ml(i,:) = crb_x * h_est'*C^-1 * y2(i,:)'; %H_ml * y2(i,:)';
                 
                 X_est_mmse(i,:) = (h_est'*C^-1*h_est + 1/sigma_X2(snr)*ones(4))^-1* h_est'*C^-1 * y2(i,:)';
 
                 MSE_x(Q_index , snr, p) = mean((X_est_ml(i,:) - x2(i,:)).^2);%,"all"
                 
                 MSE2_x(Q_index , snr, p) = mean((X_est_mmse(i,:) - x2(i,:)).^2);
             end
        end%MC
    end
end

MSE_h = mean(MSE_h, 3);
MSE_x = mean(MSE_x, 3);
MSE2_x = mean(MSE2_x, 3);

%channel
figure;grid on; box on;
loglog(sigma_X2, MSE_h)
%hold on;
title('MSE_h vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(floor(Q)));
title(lgd,'Q')

%X
data = (P*ones(1, length(Q)) - Q);
figure;
loglog(sigma_X2, MSE_x)
hold on;
xline(1)
title('MSE_x(ML) vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')


%X MMSE
figure;
loglog(sigma_X2, MSE2_x)
hold on;
xline(1)
title('MSE_x(MMSE) vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')


%metric

for Q_index = 1:length(Q)
    for k = 1:length(sigma_X2)
        metric(Q_index, k) = (1-Q(Q_index)/P)*log(1+sigma_X2(k)/MSE_x(Q_index, k));
    end
end

% figure;
% loglog(sigma_X2./MSE_x, metric)
% hold on;
% xline(1)
% title('Metric vs \sigma^2_x')
% xlabel('\sigma_X^2')
% ylabel('Metric')
% lgd = legend(string(Q));
% title(lgd,'Q')


figure;
plot(Q, metric')
hold on;
xline(1)
title('Metric vs Q')
xlabel('Q')
ylabel('Metric')
legend(string(SNR));
%  figure;
%  loglog(sigma_X2, MSE2_x)
%  title('MSE_x vs \sigma^2_x')
%  xlabel('\sigma_X^2')
%  ylabel('MSE_x')
%  legend(string(floor(Q)));

