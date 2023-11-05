%% %% 1.3 MIMO deconvolution
clear all; %clc;
M = 4;
N = 4;

SNR = -10:2:30; %dB
P = 100;
sigma_X2 = 10.^(SNR/10);

Q = floor(linspace(4, P-1, 20));%Q_min>=K(M-1)+1

rho = 0.1;
alpha = 0.5;
mse_h_ml = zeros(length(sigma_X2), length(Q));
mse_x_ml = zeros(length(sigma_X2), length(Q));
mse_x_mmse = zeros(length(sigma_X2), length(Q));

MC = 100;
h = zeros(M, N);
for l= 1:N
    for i = 1:M
        h(i,l) = alpha^abs(i-l);
    end
end

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/199'])
    for snr = 1:length(sigma_X2)
        for p =1 :MC%MC
            x =  randn(floor(Q(Q_index)),N);%sigma_X2(snr).*
            c = Cov_True(rho, M)/ sigma_X2(snr);%4X4
            w = chol(c, 'lower')* randn(M, floor(Q(Q_index)));%%%%%%%%%%%%%%%%%%%%

            W = (reshape(w',1,[]))';
            H = (reshape(h,1,[]))';
            X = kron(eye(M), x);
            Y = X*H+W;
            %Y = X*H;
                
%             r = zeros(1, M*Q(Q_index));
%             r(1) = 1;
%             for u = 1:length(r)
%                 if mod(u, Q(Q_index)+1) == 0
%                     r(u) = rho;
%                 end
%             end
% 
%             C = toeplitz(r)/ sigma_X2(snr);
            C = kron(eye(floor(Q(Q_index))), c);
            %crb = (X'*X)^-1;
            %H_EST = crb*X'*Y;
            crb = (X'*C^-1*X)^-1;
            H_EST = crb*X'*C^-1*Y;
            mse_h_ml(snr,Q_index) = mse_h_ml(snr,Q_index) + mean((h(:) - H_EST).^2)/MC;

            %MSE_h(Q_index , snr, p) = mean((H_EST - H).^2,"all");
            

            x2 =  randn(N, floor(P-Q(Q_index)));%sigma_X2(snr).*
            w2 = chol(c, 'lower')* randn(M, floor(P-Q(Q_index)));
            y2 = h * x2 + w2;
            H_ml = reshape(H_EST, M, N);
            X_EST_ML = (H_ml' * c^-1 * H_ml) \ H_ml' * c^-1 * y2;
            X_EST_MMSE = (H_ml' * c^-1 * H_ml + eye(N)) \ H_ml' * c^-1 * y2;
%             X2 =  reshape(x2,[],1);
%             h2 = reshape(H_EST,4,[]);
%             W2 =  reshape(w2,[],1);
            mse_x_ml(snr,Q_index) = mse_x_ml(snr,Q_index) + mean((x2(:) - X_EST_ML(:)).^2)/MC;
            mse_x_mmse(snr,Q_index) = mse_x_mmse(snr,Q_index) + mean((x2(:) - X_EST_MMSE(:)).^2)/MC;
% H2 =kron(eye(P-Q(Q_index)), h2);
%             Y2 = H2*X2+W2;
%             %Y2 = H2*X2;
%               
% %             r2 = zeros(1, M*(P-Q(Q_index)));
% %             r2(1) = 1;r2(2) = rho; r2(3) = rho; r2(4) = rho;
% %             C2 = toeplitz(r2);
%             crb_x = (H2'*C2^-1 * H2)^-1;
%             X_EST_ML = crb_x*H2'*C2^-1*Y2;
%             %crb_x = (H2'* H2)^-1;
%             %X_EST_ML = crb_x*H2'*Y2;
%             MSE_x(Q_index , snr, p) = mean((X_EST_ML - X2).^2,"all");%,"all"
%             %crb_x = (h_est' * h_est)^-1;%crb_x = (h_est' * C^-1 * h_est)^-1;
            %X_est_ml = crb_x * h_est'*C^-1* y2; % X_est_ml(i,:) = crb_x * h_est'*C^-1 * y2(i,:)';
            %X_est_ml = crb_x * h_est'*  y2;
            % Xin = (h_est'*h_est + (1/sigma_X2(snr))*eye(4))^-1;
            % X_est_mmse(i,:) = Xin * h_est' * y2(i,:)';% X_est_mmse(i,:) = (h_est'*C^-1*h_est + 1/sigma_X2(snr)*ones(4))^-1* h_est'*C^-1 * y2(i,:)';

           
             
            % MSE2_x(Q_index , snr, p) = mean((X_est_mmse(i,:) - x2(i,:)).^2)
        end%MC
    end
end

% MSE_h = mean(MSE_h, 3);
% MSE_x = mean(MSE_x, 3);
% %MSE2_x = mean(MSE2_x, 3);
% 
% figure;grid on; box on;
% loglog(sigma_X2, MSE_h)
% 
% 
% figure;
% loglog(sigma_X2, MSE_x)
figure;
loglog(sigma_X2, mse_h_ml);
title('MSE_h vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(floor(Q)));
title(lgd,'Q')


figure;
loglog(sigma_X2, mse_x_ml);
title('MSE_x(ML) vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')

figure;
loglog(sigma_X2, mse_x_mmse);
title('MSE_x(MMSE) vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')


for Q_index = 1:length(Q)
    for k = 1:length(sigma_X2)
        metric_ml(Q_index, k) = (1-Q(Q_index)/P)*log(1+sigma_X2(k)./mse_x_ml(k, Q_index));
        metric_mmse(Q_index, k) = (1-Q(Q_index)/P)*log(1+sigma_X2(k)./mse_x_mmse(k, Q_index));
    end
end
%figure; 
%loglog(sigma_X2, mse_x_mmse)
figure;
plot(Q, metric_ml)
title('Metric(ML) vs Q')
xlabel('Q')
ylabel('Metric')
lgd = legend(string(SNR));
title(lgd,'\sigma_X^2')

figure;
plot(Q, metric_mmse)
title('Metric(MMSE) vs Q')
xlabel('Q')
ylabel('Metric')
lgd = legend(string(SNR));
title(lgd,'\sigma_X^2')