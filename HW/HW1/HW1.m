%homework1 MIMO System Identificatiom=n and Deconvolution
%1.1Noise Generation
clear all; clc;
%defined metric SSE (SUM OF SQUARE ERRORS)
K = floor(linspace(10, 200,10));%numbere of samples
rho = linspace(0, 0.99, 10);%

M = 4;
Nmc = 1000;
% Noise parameters:
C = zeros(M,M);
SSE = zeros(length(K), length(rho),Nmc);% Metric adopted - SUM OF SQUARED ERRORS

for index_K = 1:length(K)  
    for index_rho =1:length(rho)%different rho values
        for n =1:Nmc % Monte-Carlo simulation
            C = Cov_True(rho(index_rho), M);%TRUE COV
            U = randn(M, K(index_K)); % White noise
            w_s = chol(C,"lower")* U;% Cholesky method
            C_s = cov(w_s.');%SAMPLE COV

            %%SSE%%%
            SSE(index_K,index_rho, n) =  sum(sum((C_s-C).^2));
        end
    end            
end

figure; 
grid on; box on;
semilogy(K, mean(SSE,3))
title('SSE vs K')
xlabel('K')
ylabel('SSE')
legend(string(rho));
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
    X = zeros(floor(Q(Q_index)),M);
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
            for p =1 :50%MC(>500)
                x = sigma_X2(snr).* randn(N, floor(Q(Q_index)));
                C = Cov_True(rho, floor(Q(Q_index)));
                w = chol(C, 'lower')* randn(floor(Q(Q_index)), M);%%%%%%%%%%%%%%%%%%%%
    
                for i=1:N %Creating X
                    X(:,i) = x(i,:);
                end
                
                y = zeros(floor(Q(Q_index)+1-1), M);
                for j=1:M
                    y(:,j) = X*h(j,:)' + w(:,j);
                    
                    crb = (X'*C^-1*X)^-1;
                    h_est(j,:) = (crb*X'*C^-1*y(:,j))';
    
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
legend(string(floor(Q)));


%Q value
Que=5;%4    15    27    38    50
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(Que,:,:),[length(alpha),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB(Que,:))
title('MSE_h vs \sigma^2_x',['Q =',num2str(floor(Q(Que)))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
legend(string(alpha));
%% %% 1.2(ii) MIMO estimation
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

rho = 0.01;

for Q_index = 1:length(Q)
    X = zeros(Q(Q_index),M);
    for b =1 :length(beta)
        for a = 1:length(alpha)%different alpha
            h = zeros(M, N);
            for l= 1:N
                for i = 1:M
                    for k = 1:K
                        h(i,l,k) = alpha(a)^abs(i-l)*beta(b)^k;
                    end
                end
            end
        
            for snr = 1:length(sigma_X2)
                for p =1 :100%MC
                    x = sigma_X2(snr).* randn(N,Q(Q_index));
                    C = Cov_True(rho, Q(Q_index)+K-1);
                    w = chol(C, 'lower')* randn(Q(Q_index)+K-1, M);
    
                    X = convmtx(x(1,:)',K);
                    for l=2:N %Creating X
                        conv_mat = convmtx(x(l,:)',K);
                        X = [X, conv_mat];
                    end
                                   
                    for j=1:M
                        hh = reshape(h(1,j,:),[1,K]);
                        for l =2:N
                            hh = [hh, reshape(h(1,j,:),[1,K])];
                        end
                        
                        y = zeros(Q(Q_index)+K-1, M);
                        y(:,j) = X*hh' + w(:,j);
                        
                         crb = (X'*C^-1*X)^-1;
                         h_est(j,:) = (crb*X'*C^-1*y(:,j))';
        
                         MSE_h(Q_index,b, a , snr, p) = mean((h_est(j,:) - hh).^2,"all");
                         
                         for i = 1:M
                             value(i) =  crb(i,i);
                         end
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
legend(string(Q));

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
legend(string(alpha));

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
legend(string(beta));
%% 1.3 MIMO deconvolution
clear all; clc;
M = 4;
N = 4;

SNR = -10:2:30; %dB
P = 200;
sigma_X2 = 10.^(SNR/10);

Q = floor(linspace(4, P-1, 20));%floor(linspace(4, P, 20));%Q_min>=K(M-1)+1

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
    X = zeros(floor(Q(Q_index)),M);
    for snr = 1:length(sigma_X2)
        for p =1 :20%MC
            x = sigma_X2(snr).* randn(N, floor(Q(Q_index)));
            C = Cov_True(rho, floor(Q(Q_index)));%4X4
            w = chol(C, 'lower')* randn(floor(Q(Q_index)), M);%%%%%%%%%%%%%%%%%%%%

            for i=1:N %Creating X
                X(:,i) = x(i,:);
            end
           
            y = zeros(floor(Q(Q_index)+1-1), M);
            for j=1:M
                y(:,j) = X*h(j,:)' + w(:,j);
                
                %X_ml = (X.' / C *X) \ X.' / C;
                crb = (X'*C^-1*X)^-1;
                h_est(j,:) = (crb*X'*C^-1*y(:,j))';%(X_ml*y(:,j))';

                MSE_h(Q_index , snr, p) = mean((h_est(j,:) - h(j,:)).^2, "all");%,"all"
            end
            
             x2 = sigma_X2(snr).* randn(floor(P-Q(Q_index)), N);
             C2 = Cov_True(rho, floor(P-Q(Q_index)));
             w2 = chol(C2, 'lower')* randn(floor(P-Q(Q_index)), M);
             
 
             y2 = zeros(floor(P-Q(Q_index)), M);
             for i = 1:floor(P-Q(Q_index))
                 y2(i,:) = (h*x2(i, :)' + w2(i,:)')';
                 
                 %%H_ml = (h_est.'*h_est) \ h_est.';%uncorrelated between each other
                 
                 %H_ml = (h_est.' /Cov_True(rho, M)*h_est) \ h_est.'* Cov_True(rho, M);
                 crb_h = (h_est' * Cov_True(rho, M)*h_est)^-1; h_est.'* Cov_True(rho, M);
                 X_est(i,:) = crb_h * h_est' * y2(i,:)'; %H_ml * y2(i,:)';
                 
                 MMSE_x(i,:) = (h_est.' / Cov_True(rho, M) \ h_est + 1/sigma_X2(snr)*ones(4))\ h_est.'/ Cov_True(rho, M)* y2(i,:)';
 
                 MSE1_x(Q_index , snr, p) = mean((X_est(i,:) - x2(i,:)).^2);%,"all"
                 MSE2_x(Q_index , snr, p) = mean((MMSE_x(i,:) - x2(i,:)).^2);
             end


        end%MC
    end
end
MSE_h = mean(MSE_h, 3);
 MSE1_x = mean(MSE1_x, 3);
 MSE2_x = mean(MSE2_x, 3);
%CRB=mean(CRB, 3);


 %alpha value
 figure;grid on; box on;
 loglog(sigma_X2, MSE_h)
 %hold on;
 title('MSE_h vs \sigma^2_x')
 xlabel('\sigma_X^2')
 ylabel('MSE_h')
 legend(string(floor(Q)));

 data = (P*ones(1, length(Q)) - Q);
 figure;
 loglog(sigma_X2, MSE1_x)
 title('MSE_x vs \sigma^2_x')
 xlabel('\sigma_X^2')
 ylabel('MSE_x')
 legend(string(floor(Q)));
 %legend(string(data));
 
 figure;
 loglog(sigma_X2, MSE2_x)
 title('MSE_x vs \sigma^2_x')
 xlabel('\sigma_X^2')
 ylabel('MSE_x')
 legend(string(floor(Q)));

 
%COMPUTING MATRIX H to achieve the linearity           
%             H1 = h_est(1,1)*eye(P-Q(Q_index));H2 = h_est(2,1)*eye(P-Q(Q_index));H3 = h_est(3,1)*eye(P-Q(Q_index));H4 = h_est(4,1)*eye(P-Q(Q_index));
%             for i = 2:N
%                 H1 = [H1, h_est(1,i)*eye((P-Q(Q_index)))];
%                 H2 = [H2, h_est(2,i)*eye((P-Q(Q_index)))];
%                 H3 = [H3, h_est(3,i)*eye((P-Q(Q_index)))];
%                 H4 = [H4, h_est(4,i)*eye((P-Q(Q_index)))];
%             end
% 
%             H1 = [H1, zeros(P-Q(Q_index), 4*(P-Q(Q_index))),zeros(P-Q(Q_index), 4*(P-Q(Q_index))),zeros(P-Q(Q_index), 4*(P-Q(Q_index)))];
%             H2 = [zeros(P-Q(Q_index), 4*(P-Q(Q_index))), H2, zeros(P-Q(Q_index), 4*(P-Q(Q_index))), zeros(P-Q(Q_index), 4*(P-Q(Q_index)))];
%             H3 = [zeros(P-Q(Q_index), 4*(P-Q(Q_index))), zeros(P-Q(Q_index), 4*(P-Q(Q_index))), H3, zeros(P-Q(Q_index), 4*(P-Q(Q_index)))];
%             H4 = [zeros(P-Q(Q_index), 4*(P-Q(Q_index))), zeros(P-Q(Q_index), 4*(P-Q(Q_index))), zeros(P-Q(Q_index), 4*(P-Q(Q_index))), H4];
%             
%             H = zeros(4*(P-Q(Q_index)), 4*4*(P-Q(Q_index)));% 4*(P-Q + K -1) X 4(P-Q)
%             H = cat(1, H1, H2, H3, H4);
%             %H 
%             
%             x2 = sigma_X2(snr).* randn(N, floor(P-Q(Q_index))+1-1);
%             c2 = Cov_True(rho, floor(P-Q(Q_index))+1-1);
%             w = chol(c2, 'lower')* randn(floor(P-Q(Q_index))+1-1, M);%%%%%%%%%%%%%%%%%%%%
%            
%             X2 = cat(1, reshape(x2, [],1), reshape(x2, [],1), reshape(x2, [],1), reshape(x2, [],1));%%%%%ERRROR
%            
%             y2 = H*X2;
%             C2 = Cov_True(rho, floor(4*(P-Q(Q_index))+1-1));
%             H_ml = (H.' /C2*H) \ H.'* C2;
%             X_est = H_ml * y2;

            
            %CRB(Q_index, snr, p) = mean(value);
            %disp(['MC iteration =', p, '/200', '   Q =', Q(Q_index), '/199'])