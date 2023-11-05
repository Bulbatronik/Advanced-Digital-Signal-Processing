%% hw2 Array Processing
clc, clear, close all;
load('HW2.1.mat');
%% 1)single source tracking
[N, L, P] = size(y);
SNR_dB = 0:2:40;
SNR = 10.^(SNR_dB/10);
fo = 2.9e9;
c = 3*10^8;
d = c / (2*fo);
dl = 50;

K = L/dl;% number of steps considered

lo = c / fo;
wo = 2*pi/lo;
V = 2;
dt = 1;

aoa_true = rad2deg(atan2(track_s1(2,:), track_s1(1,:)));
%% See the trajectory
pos_a = [(0:N-1).'*d zeros(N,1) ];

figure; hold on;
plot(track_s1(1,:), track_s1(2,:), '-x');
plot(pos_a(:,1), pos_a(:,2), 'o');
%% (causes problems)

A = [1 0 dt 0; ...
     0 1 0 dt; ...
     0 0 1 0; ...
     0 0 0 1];

C_alfa = [0 0 0 0; ...
          0 0 0 0; ...
          0 0 2*cos(pi/360) 0; ...
          0 0 0 2*sin(pi/360)] + abs(1e-3*diag(randn(4,1)));

B=@(x)[1/(1+(x(2)/(x(1)))^2) *(-x(2))/((x(1)^2)),...
       1/(1+(x(2)/(x(1)))^2) *(1/(x(1))),...
       0,...
       0];

az = 0: pi/360 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);

angle_pred = zeros(K,P);

beam_ang = zeros(1,K);
for p = 1:P
    angle_pred(1,p) = rad2deg(atan2(track_s1(2), track_s1(1)));

    %initialization
    theta_pp = ([track_s1(1), track_s1(2), track_s1(3)*cos(track_s1(4))*dt, track_s1(3)*sin(track_s1(4))*dt])';
    P_pp = eye(4)* 10^3;
    
    %state_est = zeros(4,K, P);
    %state_est(:,1,p) = theta_pp;
    measurement(1, p) = atan2(track_s1(2), track_s1(1));

    for i = 2:K
        %Find the angle/observation
        samples = 1+(i-1)*dl:i*dl;
        %angular_spect = mean(abs(b_c' * y(:,samples, p)).^2, 2);
        
        [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
        %beam_ang(i) = az(ind);
    
        %Prediction
        theta_cp = A * theta_pp;
        P_cp = A * P_pp * A' + C_alfa;
        
        % Update:
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
        b = atan2(theta_cp(2), theta_cp(1));

        angle_pred(i,p) = rad2deg(b);

        theta_cc = theta_cp + G*(az(ind)-b);
        measurement(i, p) = az(ind);
        P_cc = P_cp - G *B(theta_cp)*P_cp;
        
        %state_est(:,i,p) = theta_cc;
    
        theta_pp = theta_cc;
        P_pp = P_cc;
        
    end
end
%% MSE(SMOOTH)
MSE =zeros(P,1);
for p = 1:P
  MSE(p) = mean( (aoa_true' - angle_pred(:, p)).^2 );
end

figure;
plot(SNR_dB, 10*log10(MSE));

%% close all
figure;
plot(aoa_true,'-');
grid on;
hold on
plot(angle_pred(:,11), '-x')
plot(angle_pred(:,21), '-o')

plot(rad2deg(measurement(:,1)), 'x')
plot(rad2deg(measurement(:,21)), 'o')
%% 
% figure;
% plot(track_s1(1,:),track_s1(2,:), '-x');
% hold on
% 
% plot(state_est(1,:),state_est(2,:), '-o');
%% A,B - non-linear (x,y,v,psi)
a = @(x)[x(1) + x(3)*cos(x(4)) * dt; ... X
         x(2) + x(3)*sin(x(4)) * dt; ... Y
         x(3);                ... V
         x(4)];    ... Psi 

A=@(x)[1, 0, cos(x(4)) * dt, -x(3)*sin(x(4))*dt;...
       0, 1, sin(x(4)) * dt, x(3)*cos(x(4))*dt;...
       0, 0, 1, 0; ...
       0, 0, 0, 1];
              
C_alfa = [0 0 0 0; ...
          0 0 0 0; ...
          0 0 0 0; ...
          0 0 0 pi/360] + abs(1e-3*diag(randn(4,1)));
 
B=@(x)[1/(1+(x(2)/(x(1)))^2) *(-x(2))/((x(1)^2)),...
       1/(1+(x(2)/(x(1)))^2) *(1/(x(1))),...
       0,...
       0];

az = 0: pi/360 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);

angle_pred = zeros(K, P);
for p = 1:P
    angle_pred(1,p) = rad2deg(atan2(track_s1(2), track_s1(1)));
    
    %initialization
    theta_pp = track_s1(:,1);
    P_pp = eye(4)* 10^3;
    
    %state_est = zeros(4,100);
    %state_est(:,1) = theta_pp;
    %save_angle(1) = rad2deg(atan2(theta_pp(2), theta_pp(1)));
    
    for i = 2 : K
        
        %Find the angle/observation
        samples = 1+(i-1)*dl:i*dl;
        [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
        %save_angle(i) = rad2deg(az(ind));
    
        %Prediction
        theta_cp = a(theta_pp);
        P_cp = A(theta_cp) * P_pp * A(theta_cp)' + C_alfa;
        
        % Update:
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
        b = atan2(theta_cp(2), theta_cp(1));

        angle_pred(i,p) = rad2deg(b);

        theta_cc = theta_cp + G*(az(ind)-b);
        P_cc = P_cp - G *B(theta_cp)*P_cp;
        
        
        theta_pp = theta_cc;
        P_pp = P_cc;
        
        %state_est(:,i) = theta_cc;
    
    end
end
%% %% MSE(SHARP)
MSE =zeros(P,1);
for p = 1:P
   MSE(p) = mean((aoa_true' - angle_pred(:, p)).^2);%measurement error add
end

figure;
plot(SNR_dB, 10*log10(MSE));
%% 
figure;
plot(aoa_true,'b-');
grid on;
hold on
plot(angle_pred(:,10), '-x')
hold on
plot(angle_pred(:,21), '-o')
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Source separation%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all;
load('HW2.2.mat');
%% 
[N, L, P] = size(y);
SNR_dB = 0:2:40;
SNR = 10.^(SNR_dB/10);
fo = 2.9e9;
c = 3*10^8;
d = c / (2*fo);
dl = 50;

K = L/dl;% number of steps considered

lo = c / fo;
wo = 2*pi/lo;
V = 2;
dt = 1;

pos_a = [(0:N-1).'*d zeros(N,1) ];

figure;
plot(track_s1(1,:),track_s1(2,:), 'r-x');
hold on
plot(track_s2(1,:),track_s2(2,:), 'b-o');
plot(pos_a(:,1), pos_a(:,2), 'o');
%% 
az = 0: pi/180 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);

p = 10;
i = 1;
samples = 1+(i-1)*dl:i*dl;
g = mean(abs(b_c' * y(:,samples, p)).^2, 2);
[value, ind] = max(g);

plot(rad2deg(az), g )



delta = 5;%*pi/180;
%az_deg = rad2deg(az);

a = (rad2deg(az(ind))-delta):1:(rad2deg(az(ind))+1+delta);
b_c = exp(1i * wo * d * (0:N-1).' * cos(deg2rad(a))) / sqrt(N);
g2 = mean(abs(b_c' * y(:,samples, p)).^2, 2);
[value, ind2] = max(g2);
figure;
plot(a, g2 )
%[value, in] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
%% 

A = [1 0 dt 0; ...
     0 1 0 dt; ...
     0 0 1 0; ...
     0 0 0 1];

C_alfa = [0 0 0 0; ...
          0 0 0 0; ...
          0 0 2*cos(pi/360) 0; ...
          0 0 0 2*sin(pi/360)];% + abs(1e-3*diag(randn(4,1)));

B=@(x)[1/(1+(x(2)/(x(1)))^2) *(-x(2))/((x(1)^2)),...
       1/(1+(x(2)/(x(1)))^2) *(1/(x(1))),...
       0,...
       0];

aoa_true = rad2deg(atan2(track_s1(2,:), track_s1(1,:)));

%delta = 2*pi/180;
az = 0: pi/180 : pi;
az_deg = rad2deg(az);

b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
delta = 5;%area around the predicted one

angle_pred = zeros(K,P);

beam_ang = zeros(1,K);
for p = 1:P
    angle_pred(1,p) = rad2deg(atan2(track_s1(2), track_s1(1)));
    
    %scan everything(but track source 1) 
    b_c_init = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
    
    start_samples = 1:50;
    g = mean(abs(b_c_init' * y(:,start_samples, p)).^2, 2);
    [value, ind] = max(g);
    x = az(ind);

    %initialization
    theta_pp = ([track_s1(1), track_s1(2), track_s1(3)*cos(track_s1(4))*dt, track_s1(3)*sin(track_s1(4))*dt])';
    P_pp = eye(4)* 10^3;
    
    %state_est = zeros(4,K, P);
    %state_est(:,1,p) = theta_pp;
    

    for i = 2:K
        %Find the angle/observation
%         samples = 1+(i-1)*dl:i*dl;
%         [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
        %beam_ang(i) = az(ind);
    
        %Prediction
        theta_cp = A * theta_pp;
        P_cp = A * P_pp * A' + C_alfa;
        
        % Update:
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
        b = atan2(theta_cp(2), theta_cp(1));

        angle_pred(i,p) = rad2deg(b);
        
        %x = az(ind);
        theta_cc = theta_cp + G*(x-b);
        P_cc = P_cp - G *B(theta_cp)*P_cp;
        
        %state_est(:,i,p) = theta_cc;
    
        theta_pp = theta_cc;
        P_pp = P_cc;
        
        %searching around the prediction
        samples = 1+(i-1)*dl:i*dl;
        point = floor(rad2deg(b));
        a = (rad2deg(point)-delta):1:(rad2deg(point)+1+delta);
        b_c = exp(1i * wo * d * (0:N-1).' * cos(deg2rad(a))) / sqrt(N);
      
        [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
        x = az(ind);
    end
end