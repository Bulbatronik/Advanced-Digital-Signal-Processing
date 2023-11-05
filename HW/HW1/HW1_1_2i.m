%% 1.2(i) MIMO estimation
clear all; clc;
M = 4;
N = 4;
%K =1; 
%b = 1;
SNR = -10:2:30; %dB
%SNR_linear 
sigma_X2 = 10.^(SNR/10);
%sigma_X2 = -10:2:30; %dB
alpha = linspace(0,.99,5);
Q = linspace(5, 50, 5);%Q_min>=K(M-1)+1
rho = 0.1;

for Q_index = 1:length(Q)
    %X = zeros(floor(Q(Q_index)),M);
    for a = 1:length(alpha)%different alpha
       % h = h_create(M, N, K, alpha(a), b);%create the h matrix(for a given alpha)
        h = zeros(M, N);
        for l= 1:N
            for i = 1:M
                %for k = 1:K
                    h(i,l) = alpha(a)^abs(i-l);%*b^k;
                %
            end
        end
    
        for snr = 1:length(sigma_X2)
            for p =1 :500%MC(>500)
                x = sigma_X2(snr).* randn(floor(Q(Q_index)), N);
                C = Cov_True(rho, M);
                w = chol(C, 'lower')* randn(M, floor(Q(Q_index)));%%%%%%%%%%%%%%%%%%%%
    
                
                y = zeros(floor(Q(Q_index)+1-1), M);
                for j=1:M
                    y(:,j) = x*h(j,:)' + w(j,:)';
                    
                    crb = (x'*x)^-1;
                    h_est(j,:) = (crb*x'*y(:,j))';
    
                    MSE_h(Q_index, a , snr, p) = mean((h_est(j,:) - h(j,:)).^2,"all");

                    value =zeros(M,1);
                    for i = 1:M
                        value(i) =  crb(i,i);
                    end
                end
                CRB(Q_index, snr, p) = mean(value);
            end%MC
        end
    end
end
MSE_h = mean(MSE_h, 4);
CRB=mean(CRB, 3);


%alpha value
alph=3;%0    0.2475    0.4950    0.7425    0.9900
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(:,alph,:),[length(Q),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB)
title('MSE_h vs \sigma^2_x',['\alpha =',num2str(alpha(alph))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(floor(Q)));
title(lgd,'Q')


%Q value
Que=4;%4    15    27    38    50
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(Que,:,:),[length(alpha),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB(Que,:))
title('MSE_h vs \sigma^2_x',['Q =',num2str(floor(Q(Que)))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(alpha));
title(lgd,'\alpha')