%EX 3 Classififcation
clc, clear, close all;
load('HW2.3.mat');
%% 
[N, L, P] = size(y1);
SNR_dB = 0:5:40;
SNR = 10.^(SNR_dB/10);

fo = 2.9e9;
c = 3*10^8;
d = c / (2*fo);
dl = 50;

%K = L/dl;% number of steps considered

lo = c / fo;
wo = 2*pi/lo;
%V = 2;
dt = 1;

pos_a1 = [-15+(0:N-1).'*d zeros(N,1) ];
pos_a2 = [15+(0:N-1).'*d zeros(N,1) ];

 figure;
 plot(track_s1(:,1),track_s1(:,2), 'r-x');
 hold on
 plot(track_s1(1,1),track_s1(1,2), 'g-x');
 grid on
 grid minor
 plot(pos_a1(:,1), pos_a1(:,2), 'o');
 plot(pos_a2(:,1), pos_a2(:,2), 'o');
%% 
%sigmoid function
sigma=@(x)1/(1+exp(-x));
%% 

%Find the angle/observation


% 
% % xtest=xtrain;
% % ytest=ytrain;
% % 
% % a = zeros(2,401);
% angle1 = zeros(K,1);
% angle2 = zeros(K,1);
% 
 az = 0: pi/(10*360) : pi;
 b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
 steps = 100;


 p = 9;


for i =1:steps           
    samples = 1+(i-1)*dl:i*dl;
    
    [value1, ind1] = max(mean(abs(b_c' * y1(:,samples, p)).^2, 2));
    angle1(i) = rad2deg(az(ind1));
    
    [value2, ind2] = max(mean(abs(b_c' * y2(:,samples, p)).^2, 2));
    angle2(i) = rad2deg(az(ind2));
end

for i = 1:70
    X_train(i,:) = [1 angle1(i), angle2(i)];
end
 
for i = 1:30
    X_test(i,:) = [1 angle1(i+70), angle2(i+70)];
end
 %Y_train = labels(1:steps/2);
 %Y_test = labels(K/2+1:end);
 


%train - 12+12 elements (both sensors), batch = 50

%% 

mu = 0.1;
%batch_size = 50;
K = 400;% #CLASSES


% x_train = X(1:2500,:);
% x_test= X(2501:5000,:);

labels = zeros(steps,1);
for i = 1 : steps
    labels(i) = Class(track_s1(i,:));
end

a = [ones(1, K); zeros(2, K)];
for i = 1:size(X_train, 1)
    y = 0; 
    for k = 1:K
        y = 0;
        if labels(i) == k
            y = 1;
        end
        loss = 0;
        for k1 = 1:K       
            z(k1) = a(:,k1)'*X_train(i,:)';
        end
        y_e = exp(z(k))/sum(exp(z));
        a(:,k) = a(:,k) + mu*(y - y_e)*X_train(i,:)';
    end
end

Pr = zeros(400,1);
for k1 = 1:K       
    z(k1) = a(:,k1)'*X_test(1,:)';
end


for i = 1: K
    result_pred(i) = exp(z(i))/sum(exp(z));
end

[v,ind] = max(Pr);
result_pred = ind;
result_true = labels(71);

