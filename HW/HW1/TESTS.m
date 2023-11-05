%% 1.2(i) MIMO estimation
clear all; clc;
M = 4;
N = 4;
SNR = -10:2:30; %dB
%SNR_linear 
sigma_X2 = 10.^(SNR/10);

alpha = linspace(0,0.99,5);
Q = floor(linspace(4, 50, 5));%Q_min>=K(M-1)+1
rho = 0.1;

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/50'])
    %X = zeros(floor(Q(Q_index)),M);
    for a = 1:length(alpha)%different alpha
        h = zeros(M, N);
        for l= 1:N
            for i = 1:M
                h(i,l) = alpha(a)^abs(i-l);
            end
        end
    
        for snr = 1:length(sigma_X2)
            for p =1 :100%MC(>500)
                x = sigma_X2(snr).* randn(Q(Q_index),N);%NxQ
                c = Cov_True(rho, M);
                w = chol(c, 'lower')* randn(M, Q(Q_index));%%%%%%%%%%%%%%%%%%%%
                
                W = (reshape(w',1,[]))';
                H = (reshape(h,1,[]))';
                X = kron(eye(M), x);
                Y = X*H+W;
                %Y = X*H;
                    
%                 r = zeros(1, M*Q(Q_index));
%                 r(1) = 1;
%                 for u = 1:length(r)
%                     if mod(u, Q(Q_index)+1) == 0
%                         r(u) = rho;
%                     end
%                 end
% 
%                 C = toeplitz(r);
                C = kron(eye(floor(Q(Q_index))), c);
                %crb = (X'*X)^-1;
                %H_EST = crb*X'*Y;
                crb = (X'*C^-1*X)^-1;
                H_EST = crb*X'*C^-1*Y;
                MSE_h(Q_index, a , snr, p) = mean((H_EST - H).^2,"all");
                
                value =zeros(M*N,1);
                for i = 1:length(value)
                    value(i) =  crb(i,i);
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
Que=2;%4    15    27    38    50
figure;grid on; box on;
loglog(sigma_X2, reshape(MSE_h(Que,:,:),[length(alpha),length(sigma_X2)]))
hold on;
loglog(sigma_X2, CRB(Que,:))
title('MSE_h vs \sigma^2_x',['Q =',num2str(floor(Q(Que)))] )
xlabel('\sigma_X^2')
ylabel('MSE_h')
lgd = legend(string(alpha));
title(lgd,'\alpha')
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
alpha = linspace(0,0.99,5);
Q = floor(linspace(13, 50, 5));%Q_min>=K(M-1)+1

rho = 0.1;

for Q_index = 1:length(Q)
    disp( ['Q =', num2str(Q(Q_index)), '/50'])
    %X = zeros(Q(Q_index),M);
    for b =1 :length(beta)
        for a = 1:length(alpha)%different alpha
            for l= 1:N
                for i = 1:M
                    for k = 1:K
                        h(i,l,k) = alpha(a)^abs(i-l)*beta(b)^k;
                    end
                end
            end
        
            for snr = 1:length(sigma_X2)
                for p =1 :100%MC
                    x = sigma_X2(snr).* randn(N, Q(Q_index));
                    C = Cov_True(rho, M);
                    w = chol(C, 'lower')* randn(M, Q(Q_index)+K-1);
    
                    X1 = convmtx(x(1, :)',K);
                    X2 = convmtx(x(2, :)',K);
                    X3 = convmtx(x(3, :)',K);
                    X4 = convmtx(x(4, :)',K);
                    XX = [X1,X2,X3,X4];
                    X = kron(eye(M), XX);
                    %W = reshape(w,[],1);
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
                    %Y = X*H

                    r = zeros(1, M*(Q(Q_index)+K-1));
                    r(1) = 1;
                    for u = 1:length(r)
                        if mod(u, Q(Q_index)+1) == 0
                            r(u) = rho;
                        end
                    end
    
                    C = toeplitz(r);
                    %crb = (X'*X)^-1;
                    %H_EST = crb*X'*Y;
                    crb = (X'*C^-1*X)^-1;
                    H_EST = crb*X'*C^-1*Y;
                    MSE_h(Q_index,b, a , snr, p) = mean((H_EST - H).^2,"all");
                    
                    value =zeros(K*M*N,1);
                    for i = 1:length(value)
                        value(i) =  crb(i,i);
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
loglog(sigma_X2, reshape(MSE_h(:,bet1, alph,:),[length(Q), length(sigma_X2)]))
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
%% %% 1.3 MIMO deconvolution
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
                y(:,j) = x*h(j,:)'+ w(j,:)';
                
                %X_ml = (X.' / C *X) \ X.' / C;
                crb_h = (x'*x)^-1;
                h_est(j,:) = (crb_h*x'*y(:,j))';%(X_ml*y(:,j))';

                MSE_h(Q_index , snr, p) = mean((h_est(j,:) - h(j,:)).^2, "all");%,"all"
            end
            
             x2 = sigma_X2(snr).* randn(N, floor(P-Q(Q_index)));
             w2 = chol(C, 'lower')* randn(M, floor(P-Q(Q_index)));
              
             y2 = zeros(M, floor(P-Q(Q_index)));
             
             y2 = h*x2;% + w2;
              
             
             %crb_x = (h_est'*C^-1 * h_est)^-1;
            
             crb_x = (h_est' * h_est)^-1;%crb_x = (h_est' * C^-1 * h_est)^-1;
             %X_est_ml = crb_x * h_est'*C^-1* y2; % X_est_ml(i,:) = crb_x * h_est'*C^-1 * y2(i,:)';
             X_est_ml = crb_x * h_est'*  y2;
            % Xin = (h_est'*h_est + (1/sigma_X2(snr))*eye(4))^-1;
            % X_est_mmse(i,:) = Xin * h_est' * y2(i,:)';% X_est_mmse(i,:) = (h_est'*C^-1*h_est + 1/sigma_X2(snr)*ones(4))^-1* h_est'*C^-1 * y2(i,:)';

             MSE_x(Q_index , snr, p) = mean((X_est_ml - x2).^2,"all");%,"all"
             
            % MSE2_x(Q_index , snr, p) = mean((X_est_mmse(i,:) - x2(i,:)).^2)
        end%MC
    end
end

MSE_h = mean(MSE_h, 3);
MSE_x = mean(MSE_x, 3);
%MSE2_x = mean(MSE2_x, 3);

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
%data = (P*ones(1, length(Q)) - Q);
figure;
loglog(sigma_X2, MSE_x)
hold on;
xline(1)
title('MSE_x(ML) vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE_x')
lgd = legend(string(Q));
title(lgd,'Q')


% %X MMSE
% figure;
% loglog(sigma_X2, MSE2_x)
% hold on;
% xline(1)
% title('MSE_x(MMSE) vs \sigma^2_x')
% xlabel('\sigma_X^2')
% ylabel('MSE_x')
% lgd = legend(string(Q));
% title(lgd,'Q')
% 
% 
% %metric
% 
for Q_index = 1:length(Q)
    for k = 1:length(sigma_X2)
        metric(Q_index, k) = (1-Q(Q_index)/P)*log(1+sigma_X2(k)/MSE_x(Q_index, k));
    end
end


for i = length(Q)
    for j = 1: length(sigma_X2)
        v =sigma_X2(j)./MSE_x(i,j)
        g(i,j) = v;
    end
end
 

figure;
loglog(sigma_X2./MSE_x, metric)
hold on;
xline(1)
title('Metric vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('Metric')
lgd = legend(string(Q));
title(lgd,'Q')


figure;
plot(Q, metric')
hold on;
title('Metric vs Q')
xlabel('Q')
ylabel('Metric')
legend(string(SNR));


 figure;
 loglog(sigma_X2, MSE2_x)
 title('MSE_x vs \sigma^2_x')
 xlabel('\sigma_X^2')
 ylabel('MSE_x')
 legend(string(floor(Q)));


