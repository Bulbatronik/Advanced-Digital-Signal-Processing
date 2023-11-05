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
 p = 9;
 az = 0: pi/(50*360) : pi;
 b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
 b_c2 = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
 steps = 100;

for i =1:steps           
    samples = 1+(i-1)*dl:i*dl;
    
    [value1, ind1] = max(mean(abs(b_c' * y1(:,samples, p)).^2, 2));
    angle1(i) = az(ind1);
    
    [value2, ind2] = max(mean(abs(b_c2' * y2(:,samples, p)).^2, 2));
    angle2(i) = az(ind2);
end
%% A,B - non-linear (x,y,v,psi)TOO HIGH, DON'T USE
%v = 2;
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
          0 0 0 pi/180];% + abs(1e-3*diag(randn(4,1)));
 

B=@(x)[1/(1+(x(2)/(x(1)+15))^2) *(-x(2))/((x(1)^2)), 1/(1+(x(2)/(x(1)+15))^2) *(1/(x(1)+15)), 0, 0;...
       1/(1+(x(2)/(x(1)-15))^2) *(-x(2))/((x(1)^2)), 1/(1+(x(2)/(x(1)-15))^2) *(1/(x(1)-15)), 0, 0];


az = 0: pi/3600 : pi;
b_c = exp(1i * wo * d * (0:N-1).' * cos(az)) / sqrt(N);
K=100;
angle_pred = zeros(2, K, P);


sig = 0.5/2*eye(2);
mu = zeros(2, 20, 20);

for i = 1:20
    for j = 1:20
        mu(:,i,j) = [-4.75+0.5*(i-1), 29.75-0.5*(j-1)];
    end
end

for p = 1:P
    %initialization
    theta_pp = [track_s1(1,1);  track_s1(1,2); 1; 1.8];
    P_pp = eye(4)* 10^3;
    
    state_est = zeros(4,100);
    state_est(:,1) = theta_pp;
    %save_angle(1) = rad2deg(atan2(theta_pp(2), theta_pp(1)));
    
    for i = 2 : K        
        %Find the angle/observation
        samples = 1+(i-1)*dl:i*dl;
        [value1, ind1] = max(mean(abs(b_c' * y1(:,samples, p)).^2, 2));
        [value2, ind2] = max(mean(abs(b_c' * y2(:,samples, p)).^2, 2));
        %save_angle(i) = rad2deg(az(ind));
        x = [az(ind1); az(ind2)];

        %Prediction
        theta_cp = a(theta_pp);
        P_cp = A(theta_cp) * P_pp * A(theta_cp)' + C_alfa;
        
        % Update:
        G = P_cp*B(theta_cp)'*(B(theta_cp)*P_cp*B(theta_cp)'+eye(2)*.1/SNR(p))^-1;
        
        b = [atan2(theta_cp(2), theta_cp(1)+15); atan2(theta_cp(2), theta_cp(1)-15)];

        %angle_pred(i,p) = rad2deg(b);

        theta_cc = theta_cp + G*(x-b);
        P_cc = P_cp - G *B(theta_cp)*P_cp;
         
        theta_pp = theta_cc;
        P_pp = P_cc;
        
        state_est(:,i) = theta_cc;    
    end
end