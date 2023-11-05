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
plot(track_s1(1,1),track_s1(2,1), 'g-x');
plot(track_s2(1,:),track_s2(2,:), 'b-o');
plot(track_s2(1,1),track_s2(2,1), 'g-o');
plot(pos_a(:,1), pos_a(:,2), 'y-o');
%% 
az = 0: pi/1080 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);

p = 21;
i = 1;
samples = 1+(i-1)*dl:i*dl;
g = mean(abs(b_c' * y(:,samples, p)).^2, 2);
[value, ind] = max(g);

figure;
plot(rad2deg(az), g )



delta = 10;

a = (rad2deg(az(ind))-delta):1:(rad2deg(az(ind))+1+delta);
b_c = exp(1i * wo * d * (0:N-1).' * cos(deg2rad(a))) / sqrt(N);
g2 = mean(abs(b_c' * y(:,samples, p)).^2, 2);
[value, ind2] = max(g2);
figure;
plot(a, g2 )
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


    %scan everything
    scan_area = 0: pi/(20*360) : pi;
    b_init = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
    samples = 1:dl;
    g  = mean(abs(b_init' * y(:,samples, p)).^2, 2);
    %plot(rad2deg(scan_area), g )
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
        %plot(rad2deg(scan_area), g )
        [value, ind] = max(g);
        
        %multiply by the b(aoe)
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
    MSE(p) = abs(mean((s2_est(:,p) - s2).^2));
end

plot(SNR_dB, 10*log10(smooth(MSE)))
%% 
figure;
plot(aoa_true,'b-');
grid on;
hold on
plot(angle_pred(:,1), '-x')
%% 


%     %x = az(ind);
%     x = scan_area(ind);
%     disp(rad2deg(x));
% 
%     
% 
% 
%     b_c = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
%     samples = 1+(i-1)*dl:i*dl;
%     g  = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
%     [value, ind] = max(g);
%     disp(ind)


%angle_pred = zeros(K,P);

%beam_ang = zeros(1,K);
% for p = 1:P
%     angle_pred(1,p) = rad2deg(atan2(track_s2(2), track_s2(1)));
%     
%     %scan everything(but track source 1) 
%     scan_area = 0: pi/180 : pi;
%     b_c = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
%     samples = 1:dl;
%     g  = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
%     [value, ind] = max(g);
% 
%     %initialization
%     theta_pp = ([track_s2(1), track_s2(2), track_s2(3)*cos(track_s2(4))*dt, track_s2(3)*sin(track_s2(4))*dt])';
%     P_pp = eye(4)* 10^3;   
%     
%     for i = 2:K
%         b_c = exp(1i * wo * d * (0:N-1).' * cos(scan_area)) / sqrt(N);
%         samples = 1+(i-1)*dl:i*dl;
%         g  = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
%         [value, ind] = max(g);
%         disp(ind)
% 
% 
% 
% 
% 
% 
% 
% 
%         
%         %Prediction
%         %disp(rad2deg(x))
%         theta_cp = A * theta_pp;
%         P_cp = A * P_pp * A' + C_alfa;
%         
%         % Update:
%         G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+1/SNR(p))^-1;
%         b = atan2(theta_cp(2), theta_cp(1));
% 
%         angle_pred(i,p) = rad2deg(b);
%         
%         %x = az(ind);
%         x = scan_area(ind);
%         disp(rad2deg(x));
% 
%         theta_cc = theta_cp + G*(x-b);
%         P_cc = P_cp - G *B(theta_cp)*P_cp;
%         
%         %state_est(:,i,p) = theta_cc;
%     
%         theta_pp = theta_cc;
%         P_pp = P_cc;
%         
%         %searching around the prediction
%         scan_area = deg2rad(floor(rad2deg(b))-delta:floor(rad2deg(b))+delta);
% %         samples = 1+(i-1)*dl:i*dl;
% %         point = floor(rad2deg(b));
% %         a = (rad2deg(point)-delta):1:(rad2deg(point)+1+delta);
% %         b_c = exp(1i * wo * d * (0:N-1).' * cos(deg2rad(a))) / sqrt(N);
% %       
% %         [value, ind] = max(mean(abs(b_c' * y(:,samples, p)).^2, 2));
% %         x = az(ind);
%     end
% end
%% 
figure;
plot(aoa_true,'b-');
grid on;
hold on
plot(abs(angle_pred(:,1)), '-g')
hold on
plot(abs(angle_pred(:,21)), '-y')
title('Angle vs Time')
xlabel('Time')
ylabel('Angle')
legend("True observation", "SNR(dB) = 0 dB", "SNR(dB) = 40 dB")
