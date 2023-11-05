%% %% 1.3 MIMO deconvolution
clear all; %clc;
M = 4;
N = 4;

MC = 50;

SNR = -10:2:30; %dB
P = 100;
sigma_X2 = 10.^(SNR/10);

Q_vec = floor(linspace(4, P-1, 10));%Q_min>=K(M-1)+1

mse_h_ml = zeros(length(sigma_X2), length(Q_vec));
mse_x_ml = zeros(length(sigma_X2), length(Q_vec));
mse_x_mmse = zeros(length(sigma_X2), length(Q_vec));

rho = 0.1;
alpha = 0.5;

h = zeros(M, N);
for l= 1:N
    for i = 1:M
        h(i,l) = alpha^abs(i-l);
    end
end

for iq = 1 : length(Q_vec)
    
    disp( ['Q = ', num2str(Q_vec(iq)), '/', num2str(P)])
    q =  floor(Q_vec(iq));
    
    for i = 1:length(sigma_X2)
        
        for p = 1 : MC
            
            % System definition
            x = randn(N, P);
            Cw = Cov_True(rho, M) / sigma_X2(i);
            w = chol(Cw, 'lower') * randn(M, P);
            
            y = h * x + w;
            
            % channel estimation:
            X = kron(x(:,1:q)', eye(M));
            y_ch = y(:,1:q);
            
            % ML estimation
            Cwi = kron(eye(q), Cw)^-1;
            X_ml = (X'*Cwi*X)^-1 * (X' * Cwi);
            h_ml = X_ml * y_ch(:);
            
            % deconvolution:
            H_ml = reshape(h_ml, M, N);
            y_data = y(:,q+1:end);
            x_data = x(:,q+1:end);
            
            % ML:
            x_ml = (H_ml' * Cw^-1 * H_ml) \ H_ml' * Cw^-1 * y_data;
            
            % MMSE:
            Cx = eye(N);
            x_mmse = (H_ml' * Cw^-1 * H_ml + Cx) \ H_ml' * Cw^-1 * y_data;
            
            % Evaluation:
            mse_h_ml(i,iq) = mse_h_ml(i,iq) + mean((h(:) - h_ml).^2)/MC;

            mse_x_ml(i,iq) = mse_x_ml(i,iq) + mean((x_data(:) - x_ml(:)).^2)/MC;
            mse_x_mmse(i,iq) = mse_x_mmse(i,iq) + mean((x_data(:) - x_mmse(:)).^2)/MC;
          

            %W = (reshape(w',1,[]))';
            %H = (reshape(h,1,[]))';
            %X = kron(eye(M), x);
            %Y = X*H+W;
            %Y = X*H;
                
%             r = zeros(1, M*Q(Q_index));
%             r(1) = 1;
%             for u = 1:length(r)
%                 if mod(u, Q(Q_index)+1) == 0
%                     r(u) = rho;
%                 end
%             end
% 
%             C = toeplitz(r);
%             %crb = (X'*X)^-1;
%             %H_EST = crb*X'*Y;
%             crb = (X'*C^-1*X)^-1;
%             H_EST = crb*X'*C^-1*Y;
%             MSE_h(Q_index , snr, p) = mean((H_EST - H).^2,"all");
%             
% 
%             x2 = sigma_X2(snr).* randn(N, floor(P-Q(Q_index)));
%             w2 = chol(c, 'lower')* randn(M, floor(P-Q(Q_index)));
%             
%             X2 =  reshape(x2,[],1);
%             h2 = reshape(H_EST,4,[]);
%             W2 =  reshape(w2,[],1);
%             
%             H2 =kron(eye(P-Q(Q_index)), h2);
%             Y2 = H2*X2+W2;
%             %Y2 = H2*X2;
%               
%             r2 = zeros(1, M*(P-Q(Q_index)));
%             r2(1) = 1;r2(2) = rho; r2(3) = rho; r2(4) = rho;
%             C2 = toeplitz(r2);
%             crb_x = (H2'*C2^-1 * H2)^-1;
%             X_EST_ML = crb_x*H2'*C2^-1*Y2;
%             %crb_x = (H2'* H2)^-1;
%             %X_EST_ML = crb_x*H2'*Y2;
%             MSE_x(Q_index , snr, p) = mean((X_EST_ML - X2).^2,"all");%,"all"
%             %crb_x = (h_est' * h_est)^-1;%crb_x = (h_est' * C^-1 * h_est)^-1;
%             %X_est_ml = crb_x * h_est'*C^-1* y2; % X_est_ml(i,:) = crb_x * h_est'*C^-1 * y2(i,:)';
%             %X_est_ml = crb_x * h_est'*  y2;
%             % Xin = (h_est'*h_est + (1/sigma_X2(snr))*eye(4))^-1;
%             % X_est_mmse(i,:) = Xin * h_est' * y2(i,:)';% X_est_mmse(i,:) = (h_est'*C^-1*h_est + 1/sigma_X2(snr)*ones(4))^-1* h_est'*C^-1 * y2(i,:)';
% 
%            
%              
%             % MSE2_x(Q_index , snr, p) = mean((X_est_mmse(i,:) - x2(i,:)).^2)
        end%MC
    end
end

%
figure; 

%plot(SNR, mse_x_ml);
%plot(SNR, mse_x_mmse);
loglog(sigma_X2, mse_x_ml);
figure; 
loglog(sigma_X2, mse_x_mmse)

figure; 
semilogy(Q_vec, mse_x_ml);
figure; 
semilogy(Q_vec, mse_x_mmse);

comm_eff_ml = (1 - Q_vec/P)  ./ mse_x_ml;
comm_eff_mmse = (1 - Q_vec/P) ./ mse_x_mmse;

figure; plot(Q_vec, comm_eff_ml);
figure; plot(Q_vec, comm_eff_mmse);

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