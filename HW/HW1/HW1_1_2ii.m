%% 1.2(ii) MIMO estimation
clear all; clc;
M = 4;
N = 4;
K =4; 
beta = [0.9, 0.5, 0.1];
SNR = -10:2:30; %dB
%SNR_linear 
sigma_X2 = 10.^(SNR/10);
%sigma_X2 = -10:2:30; %dB
alpha = linspace(0.1,.99,5);
Q = floor(linspace(13, 50, 5));%Q_min>=K(M-1)+1

rho = 0.1;

for Q_index = 1:length(Q)
    X = zeros(Q(Q_index),M);
    for b =1 :length(beta)
        for a = 1:length(alpha)%different alpha
            %h = zeros(M, N);
            for l= 1:N
                for i = 1:M
                    for k = 1:K
                        h(i,l,k) = alpha(a)^abs(i-l)*beta(b)^k;
                    end
                end
            end
        
            for snr = 1:length(sigma_X2)
                for p =1 :50%MC
                    x = sigma_X2(snr).* randn(Q(Q_index),N);
                    C = Cov_True(rho, M);
                    w = chol(C, 'lower')* randn(M, Q(Q_index)+K-1);
    
                    X = convmtx(x(:, 1),K);
                    for l=2:N %Creating X
                        conv_mat = convmtx(x(:,l),K);
                        X = [X, conv_mat];
                    end
                                   
                    for j=1:M
%                         for i = 1:N
                            hh = reshape(h(j,1,:),[1,K]);
                            for l =2:N
                                hh = [hh, reshape(h(j,l,:),[1,K])];
                            end
                            
                            y = zeros(Q(Q_index)+K-1, M);
                            y(:,j) = X*hh' + w(j,:)';
                            
                             crb = (X'*X)^-1;
                             h_est(:,j) = (crb*X'*y(:,j))';
            
                             MSE_h(Q_index,b, a , snr, p) = mean((h_est(:,j) - hh').^2,"all");
                             
                             for s = 1:M
                                 value(s) =  crb(s,s);
                             end
%                         end
                    end
                    CRB(Q_index, snr, p) = mean(value);
                end%MC
            end
        end
    end
end

MSE_h = mean(MSE_h, 5);
CRB=mean(CRB, 3);

%alpha value
alph=3;%0    0.2475    0.4950    0.7425    0.9900
bet1 = 1;%[0.9, 0.5, 0.1]
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(:,bet1, alph,:),[length(Q),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB)
title('MSE_h vs \sigma^2_x',['\alpha =',num2str(alpha(alph)), '\beta =', num2str(beta(bet1))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(Q));
title(lgd,'Q')

%Q value
Que=3;%4    15    27    38    50
bet2 = 1;%[0.9, 0.5, 0.1]
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(Que,bet2,:,:),[length(alpha),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB(Que,:))
title('MSE_h vs \sigma^2_x',['Q =',num2str(Q(Que)), '\beta =',num2str(beta(bet2))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(alpha));
title(lgd,'\alpha')

%beta value
Que=4;%4    15    27    38    50
alph=3;%0    0.2475    0.4950    0.7425    0.9900
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(Que,:,alph,:),[length(beta),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB(Que,:))
title('MSE_h vs \sigma^2_x',['Q =',num2str(Q(Que)), '\alpha =',num2str(alpha(alph))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(beta));
title(lgd,'\beta')