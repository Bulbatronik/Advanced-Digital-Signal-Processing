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
plot(track_s1(1,1), track_s1(2,1), 'g-o');
plot(pos_a(:,1), pos_a(:,2), 'o');
%% EKF

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

az = 0: pi/1080 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);

angle_pred = zeros(K,P);

beam_ang = zeros(1,K);
state_est = zeros(4,K, P);
for p = 1:P
    angle_pred(1,p) = rad2deg(atan2(track_s1(2), track_s1(1)));

    %initialization
    theta_pp = ([track_s1(1), track_s1(2), track_s1(3)*cos(track_s1(4))*dt, track_s1(3)*sin(track_s1(4))*dt])';
    P_pp = eye(4)* 10^3;

    for i = 2:K
        %Find the angle/observation
        samples = 1+(i-1)*dl:i*dl;
        [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
    
        %Prediction
        theta_cp = A * theta_pp;
        P_cp = A * P_pp * A' + C_alfa;
        
        % Update:
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
        b = atan2(theta_cp(2), theta_cp(1));

        angle_pred(i,p) = rad2deg(b);

        theta_cc = theta_cp + G*(az(ind)-b);
        P_cc = P_cp - G *B(theta_cp)*P_cp;
    
        theta_pp = theta_cc;
        P_pp = P_cc;        
    end
end
%% MSE
MSE =zeros(P,1);
for p = 1:P
  MSE(p) = mean( (aoa_true' - angle_pred(:, p)).^2 );
end

figure;
plot(SNR_dB, 10*log10(MSE));
title('MSE vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE')
%% close all
figure;
plot(aoa_true,'-r');
grid on;
hold on
plot(angle_pred(:,1), '-g')
hold on
grid on
plot(angle_pred(:,21), '-b')
title('Angle vs Time')
xlabel('Time')
ylabel('Angle')
legend("True observation", "SNR(dB) = 0 dB", "SNR(dB) = 40 dB")
%% 1.2 Source separation
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
plot(track_s1(1,1),track_s1(2,1), 'g-x');
plot(track_s2(1,:),track_s2(2,:), 'b-o');
plot(track_s2(1,1),track_s2(2,1), 'g-o');
plot(pos_a(:,1), pos_a(:,2), 'y-o');
%% 
A = [1 0 dt 0; ...
     0 1 0 dt; ...
     0 0 1 0; ...
     0 0 0 1];

C_alfa = [0 0 0 0; ...
          0 0 0 0; ...
          0 0 2*cos(pi/360) 0; ...
          0 0 0 2*sin(pi/360)];%+ abs(1e-3*diag(randn(4,1)));% + abs(1e-3*diag(randn(4,1)));

B=@(x)[1/(1+(x(2)/(x(1)))^2) *(-x(2))/((x(1)^2)),...
       1/(1+(x(2)/(x(1)))^2) *(1/(x(1))),...
       0,...
       0];


%p = 21;
aoa_true = rad2deg(atan2(track_s2(2,:), track_s2(1,:)))';

s2_est = zeros(5000,P);


angle_pred = zeros(K,P);


for p = 1:P

    angle_pred(1,p) = rad2deg(atan2(track_s2(2), track_s2(1)));

    %scan everything for the first time
    scan_area = 0: pi/1080 : pi;
    b_init = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
    samples = 1:dl;
    g  = mean(abs(b_init' * y(:,samples, p)).^2, 2);
    [value, ind] = max(g);
    
    b2_init = exp(1i * wo * d * (0:N-1).' * cos(scan_area(ind))) / sqrt(N);
    s2_est(1:50,p) = (b2_init' * y(:,samples, p)).';
    
    
    %initialization
    theta_pp = ([track_s2(1), track_s2(2), track_s2(3)*cos(track_s2(4))*dt, track_s2(3)*sin(track_s2(4))*dt])';
    P_pp = eye(4)* 10^3;  
    
    for i = 2:K
        
        %Prediction
        theta_cp = A * theta_pp;
        P_cp = A * P_pp * A' + C_alfa;
       
        b = atan2(theta_cp(2), theta_cp(1));%predicted angle
        
        angle_pred(i,p) = rad2deg(b);

        
        %scan around the predicted angle
        delta = 5*pi/180;
        scan_area = b - delta: pi/(20*360) : b + delta;
        
        %max in the area
        b_c = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
        samples = 1+(i-1)*dl:i*dl;
        g  = mean(abs(b_c' * y(:,samples, p)).^2, 2);
        [value, ind] = max(g);
        
        %find the [b(AoA)*y']
        b2 = exp(1i * wo * d * (0:N-1).' * cos(scan_area(ind))) / sqrt(N);
        s2_est(samples,p) = (b2' * y(:,samples, p)).';
        
        %obeservation
        x = scan_area(ind);

        % Update
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
        theta_cc = theta_cp + G*(x-b);
        P_cc = P_cp - G *B(theta_cp)*P_cp;

        theta_pp = theta_cc;
        P_pp = P_cc;
    end
end
%% MSE
MSE = zeros(P,1);
for p =1:P
    MSE(p) = mean((s2_est(:,p) - s2).^2);
end

figure
plot(SNR_dB, 10*log10(smooth(MSE)))
title('MSE vs \sigma^2_x')
xlabel('\sigma_X^2')
ylabel('MSE')

