%homework1 MIMO System Identificatiom=n and Deconvolution
%1.1Noise Generation
clear; clc;
disp('===1.1===')
%defined metric SSE (SUM OF SQUARE ERRORS)
K = floor(linspace(10, 200,10));%numbere of samples
rho = linspace(0, 0.99, 10);
M = 4;
Nmc = 2*10^4;
% Noise parameters:
C = zeros(M,M);
SSE = zeros(length(K), length(rho), Nmc);% Metric adopted - SUM OF SQUARED ERRORS

for index_K = 1:length(K)  
    for index_rho =1:length(rho)%different rho values
        for n =1:Nmc % Monte-Carlo simulation
            C = Cov_True(rho(index_rho), M);%TRUE COV
            U = randn(M, K(index_K)); % White noise
            w_s = chol(C,"lower")* U;% Cholesky method
            C_s = cov(w_s.');%SAMPLE COV
            %%SSE%%%
            SSE(index_K,index_rho) =  SSE(index_K,index_rho) + sum(sum((C_s-C).^2))/Nmc;
        end
    end            
end

figure; 
grid on; box on;
semilogy(K, mean(SSE,3))
title('1.1 SSE vs K')
xlabel('K')
ylabel('SSE')
lgd = legend(string(rho));
title(lgd,'\rho')
%% 1.2(i) MIMO estimation
%clear all; clc;
disp('===1.2(i)===')
M = 4;
N = 4;
SNR_dB = -10:2:30; %dB
SNR = 10.^(SNR_dB/10);%SNR_linear 
alpha = linspace(0,0.99,5);
Q = floor(linspace(5, 50, 5));%Q_min>=K(M-1)+1
rho = 0.1;
Nmc = 1000;%2*10^3;
MSE_h= zeros(length(Q), length(alpha) , length(SNR), Nmc);
CRB = zeros (length(Q), length(SNR), Nmc);

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/50'])
    for a = 1:length(alpha)%different alpha
        h = zeros(M, N);
        for l= 1:N
            for i = 1:M
                h(i,l) = alpha(a)^abs(i-l);
            end
        end
    
        for s = 1:length(SNR)
            for p =1 :Nmc%MC
                x = randn(Q(Q_index),N);
                c = Cov_True(rho, M)/SNR(s);
                w = chol(c, 'lower')* randn(M, Q(Q_index));
                
                W = (reshape(w',1,[]))';
                H = (reshape(h,1,[]))';
                X = kron(eye(M), x);
                Y = X*H+W;
                
                r = zeros(1, M*(Q(Q_index)));
                r(1) = 1;
                for u = 1:length(r)
                    if mod(u, Q(Q_index)+1) == 0
                        r(u) = rho;
                    end
                end
    
                C = toeplitz(r)/SNR(s);
                crb = (X'*C^-1*X)^-1;
                H_EST = crb*X'*C^-1*Y;
                MSE_h(Q_index, a , s, p) = mean((H_EST - H).^2,"all");
                
                value =zeros(M*N,1);
                for i = 1:length(value)
                    value(i) =  crb(i,i);
                end
                CRB(Q_index, s, p) = mean(value);                
           end%MC
        end
    end
end
MSE_h = mean(MSE_h, 4);
CRB=mean(CRB, 3);

%alpha value
alph=3;%0    0.2475    0.4950    0.7425    0.9900
figure;grid on; box on;
%loglog(SNR, reshape(MSE_h(:,alph,:),[length(Q),length(SNR)]))
plot(SNR_dB, 10*log10(reshape(MSE_h(:,alph,:),[length(Q),length(SNR)])))
hold on;
%loglog(SNR, CRB,'--')
plot(SNR_dB, 10*log10(CRB),'--')
title('1.2(i)MSE_h vs SNR',['\alpha =',num2str(alpha(alph))] )
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(floor(Q)));
title(lgd,'Q')


%Q value
Que=3;%4    15    27    38    50
figure;grid on; box on;
%loglog(SNR, reshape(MSE_h(Que,:,:),[length(alpha),length(SNR)]))
plot(SNR_dB, 10*log10(reshape(MSE_h(Que,:,:),[length(alpha),length(SNR)])))
hold on;
%loglog(SNR, CRB(Que,:),'--')
plot(SNR_dB, 10*log10(CRB(Que,:)), '--')
title('1.2(i)MSE_h vs SNR',['Q =',num2str(floor(Q(Que)))] )
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(alpha));
title(lgd,'\alpha')
%% 1.2(ii) MIMO estimation
%clear all; clc;
disp('===1.2(ii)===')
M = 4;
N = 4;
K =4; 
beta = [0.9, 0.5, 0.1];
SNR_dB = -10:2:30; %dB
SNR = 10.^(SNR_dB/10);%SNR_linear 
alpha = linspace(0,0.99,5);
Q = floor(linspace(14, 50, 5));%Q_min>=K(M-1)+1
rho = 0.1;
Nmc = 1000;
MSE_h = zeros(length(Q), length(beta), length(alpha),length(SNR), Nmc);
CRB = zeros (length(Q), length(SNR), Nmc);

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/50'])
    for b =1 :length(beta)
        for a = 1:length(alpha)%different alpha
            for l= 1:N
                for i = 1:M
                    for k = 1:K
                        h(i,l,k) = alpha(a)^abs(i-l)*beta(b)^k;
                    end
                end
            end
        
            for s = 1:length(SNR)
                for p =1 :Nmc%MC
                    x = randn(N, Q(Q_index));
                    C = Cov_True(rho, M)/SNR(s);
                    w = chol(C, 'lower')* randn(M, Q(Q_index)+K-1);
    
                    X1 = convmtx(x(1, :)',K);
                    X2 = convmtx(x(2, :)',K);
                    X3 = convmtx(x(3, :)',K);
                    X4 = convmtx(x(4, :)',K);
                    XX = [X1,X2,X3,X4];
                    X = kron(eye(M), XX);
                    W = (reshape(w',1,[]))';
                        
                    H = zeros(K*N*M,1);
                    f = 1;
                    for i=1:M
                        for l =1:N
                            for k =1:K
                                H(f) = h(i,l,k);
                                f=f+1;
                            end
                        end
                    end

                    Y = X*H+W;
                    
                    r = zeros(1, M*(Q(Q_index)+K-1));
                    r(1) = 1;
                    for u = 1:length(r)
                        if mod(u, Q(Q_index)+1) == 0
                            r(u) = rho;
                        end
                    end
    
                    C = toeplitz(r)/SNR(s);
                    crb = (X'*C^-1*X)^-1;
                    H_EST = crb*X'*C^-1*Y;
                    MSE_h(Q_index, b, a , s, p) = mean((H_EST - H).^2,"all");
                    
                    value =zeros(K*M*N,1);
                    for i = 1:length(value)
                        value(i) =  crb(i,i);
                    end
                    CRB(Q_index, s, p) = mean(value);
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
%loglog(SNR, reshape(MSE_h(:,bet1, alph,:),[length(Q), length(SNR)]))
plot(SNR_dB, 10*log10(reshape(MSE_h(:,bet1, alph,:),[length(Q), length(SNR)])))
hold on;
%loglog(SNR, CRB,'--')
plot(SNR_dB, 10*log10(CRB), '--')
title('1.2(ii) MSE_h vs SNR',['\alpha =',num2str(alpha(alph)), '\beta =', num2str(beta(bet1))] )
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(Q));
title(lgd,'Q')

%Q value
Que=3;%4    15    27    38    50
bet2 = 1;%[0.9, 0.5, 0.1]
figure;grid on; box on;
%loglog(SNR, reshape(MSE_h(Que,bet2,:,:), [length(alpha),length(SNR)]))
plot(SNR_dB, 10*log10(reshape(MSE_h(Que,bet2,:,:), [length(alpha),length(SNR)])))
hold on;
%loglog(SNR, CRB(Que,:),'--')
plot(SNR_dB, 10*log10(CRB(Que,:)), '--')
title('1.2(ii) MSE_h vs SNR',['Q =',num2str(Q(Que)), '\beta =',num2str(beta(bet2))] )
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(alpha));
title(lgd,'\alpha')

%beta value
Que=4;%4    15    27    38    50
alph=3;%0    0.2475    0.4950    0.7425    0.9900
figure;grid on; box on;
%loglog(SNR, reshape(MSE_h(Que,:,alph,:),[length(beta), length(SNR)]))
plot(SNR_dB, 10*log10(reshape(MSE_h(Que,:,alph,:),[length(beta), length(SNR)])))
hold on;
%loglog(SNR, CRB(Que,:),'--')
plot(SNR_dB, 10*log10(CRB(Que,:)),'--')
title('1.2(ii) MSE_h vs SNR',['Q =',num2str(Q(Que)), '\alpha =',num2str(alpha(alph))] )
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(beta));
title(lgd,'\beta')
%% 1.3 MIMO deconvolution
%clear all; clc;
disp('===1.3===')
M = 4;
N = 4;
SNR_dB = -10:2:30; %dB
P = 200;
SNR = 10.^(SNR_dB/10);%SNR_linear 
Q = floor(linspace(5, P-1, 20));%Q_min>=K(M-1)+1
rho = 0.1;
alpha = 0.5;
Nmc = 1000;%10^3;
h = zeros(M, N);
MSE_h= zeros(length(SNR), length(Q), Nmc);
MSE_x_ML= zeros(length(SNR), length(Q), Nmc);
MSE_x_MMSE= zeros(length(SNR), length(Q), Nmc);

for l= 1:N
    for i = 1:M
        h(i,l) = alpha^abs(i-l);
    end
end

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/199'])
    for s = 1:length(SNR)
        for p =1 :Nmc%MC
            x =  randn(floor(Q(Q_index)),N);
            c = Cov_True(rho, M)/ SNR(s);%4X4
            w = chol(c, 'lower')* randn(M, floor(Q(Q_index)));

            W = (reshape(w',1,[]))';
            H = (reshape(h,1,[]))';
            X = kron(eye(M), x);
            Y = X*H+W;

            %C = kron(eye(floor(Q(Q_index))), c);
            r = zeros(1, M*(Q(Q_index)));
            r(1) = 1;
            for u = 1:length(r)
                if mod(u, Q(Q_index)+1) == 0
                    r(u) = rho;
                end
            end

            C = toeplitz(r)/SNR(s);
            crb = (X'*C^-1*X)^-1;
            
            H_EST = (X'*C^-1*X)^-1*X'*C^-1*Y;
            MSE_h(s, Q_index, p) = mean((H_EST - H).^2,"all");
            

            x2 =  randn(N, floor(P-Q(Q_index)));
            w2 = chol(c, 'lower')* randn(M, floor(P-Q(Q_index)));
            y2 = h * x2 + w2;
            H_ml = reshape(H_EST, M, N);
            X_EST_ML = (H_ml' * c^-1 * H_ml) \ H_ml' * c^-1 * y2;
            X_EST_MMSE = (H_ml' * c^-1 * H_ml + eye(N)) \ H_ml' * c^-1 * y2;
         
            MSE_x_ML(s,Q_index, p) = mean((x2(:) - X_EST_ML(:)).^2);
            MSE_x_MMSE(s,Q_index, p) = mean((x2(:) - X_EST_MMSE(:)).^2);
       end%MC
    end
end

 MSE_h = mean(MSE_h, 3);
 MSE_x_ML= mean(MSE_x_ML, 3);
 MSE_x_MMSE= mean(MSE_x_MMSE, 3);


figure;
%loglog(SNR, MSE_h);
plot(SNR_dB, 10*log10(MSE_h))
title('1.3 MSE_h vs SNR')
xlabel('SNR')
ylabel('MSE_h')
lgd = legend(string(floor(Q)));
title(lgd,'Q')


figure;
%loglog(SNR, MSE_x_ML);
plot(SNR_dB, 10*log10(MSE_x_ML))
title('1.3 MSE_x(ML) vs SNR')
xlabel('SNR')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')

figure;
%loglog(SNR, MSE_x_MMSE);
plot(SNR_dB, 10*log10(MSE_x_MMSE))
title('1.3 MSE_x(MMSE) vs SNR')
xlabel('SNR')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')

metric_ml = zeros(length(Q), length(SNR));
metric_mmse = zeros(length(Q), length(SNR));
for Q_index = 1:length(Q)
    for s = 1:length(SNR)
        metric_ml(Q_index, s) = (1-Q(Q_index)/P)*log(1+(SNR(s)./MSE_x_ML(s, Q_index)));
        metric_mmse(Q_index, s) = (1-Q(Q_index)/P)*log(1+(SNR(s)./MSE_x_MMSE(s, Q_index)));
    end
end

figure;
plot(Q, metric_ml)
title('1.3 Metric vs Q')
xlabel('Q')
ylabel('Metric')
lgd = legend(string(SNR_dB));
title(lgd,'SNR')

% figure;
% plot(Q, metric_mmse)
% title('1.3 Metric(MMSE) vs Q')
% xlabel('Q')
% ylabel('Metric')
% lgd = legend(string(SNR_dB));
% title(lgd,'SNR')

% xi = linspace(min(Q), max(Q), 150);
% yi = interp1(Q, metric_ml, xi, 'spline', 'extrap');
% figure; plot(Q, yi)