%% WS#3 - Solutions
clc, clear, close all;

load('WSnew.mat');

Fs = 1e6;       % 1 [MHz]
c = physconst('lightspeed');

BSpos = [0 0 10;        ...
         0 200 10;      ...
         200 200 10;    ...
         200 0 10];
     
[Ntest, Nbs, Nsig] = size(rt);

%% Toa Estimation:

Toa = zeros(Ntest, Nbs);

USFact = 1;
Fs_up = Fs * USFact;
t_true = 0 : 1/Fs : (Nsig-1)/Fs;
t_interp = 0 : 1 / Fs_up : (Nsig-1)/Fs_up;

st_up = interp1(t_true, st, t_interp);
st_up(end+USFact:end) = [];

Nmax = length(st_up);

for np = 1 : Ntest
    
    for nb = 1 : Nbs
        r_p_b = squeeze(rt(np, nb, :));
        r_up = interp1(t_true, r_p_b, t_interp);
        r_up(end+USFact:end) = [];
        r_rt_st = xcorr(r_up, st_up);
        
        [Rm, Im] = max(r_rt_st);
        Rm_m1 = r_rt_st(Im-1);
        Rm_p1 = r_rt_st(Im+1);
        T = (Rm_m1 - Rm_p1) / (2 * (Rm_m1 + Rm_p1 - 2*Rm));
        
        Toa(np, nb) = ((Im - Nmax) + T) / Fs_up;
        
    end
end

%% Position Estimation:

A = [-2*BSpos(:,1) -2*BSpos(:,2) ones(Nbs,1)];
Pest_lls = zeros(Ntest, 2);

for np = 1 : Ntest
   
    b = ((Toa(np,:) * c).^2).' - BSpos(:,1).^2 - BSpos(:,2).^2;
    x_lls = A \ b;
    Pest_lls(np,:) = x_lls(1:2);
        
end

clc;
%% Kalman Filtering
% f(x,y,vx,vy,ax,ay) = x + vx*dt + ax * dt^2/2; vx + ax * dt; ax
%                    = y + vy*dt + ay * dt^2/2; vy + ay * dt; ay

close all
T = 1;             % measurments update 1 s;
Xs = 6;
% Dynamic model:
F = [1 0 T 0 T^2/2 0; ...
     0 1 0 T 0 T^2/2; ...
     0 0 1 0 T 0;      ...
     0 0 0 1 0 T;      ...
     0 0 0 0 1 0;       ...
     0 0 0 0 0 1];
ax = 5e1 * 8.4; ay = 4e-5 * 6; 
vx = 1e-4 * 2.4; vy = 5e1 * 1.9; 
px = 1e-5 * 1.2; py = 1e5 * 1;

Q = 1e-8 * diag([px py vx vy ax ay]);

% Measurements model:
H = [1 0 0 0 0 0; ...
     0 1 0 0 0 0];
 
R = 1*diag([2 2]);

% Initialization:
Xk = [Pest_lls(1,:) -.001 .02 0.0002 0.001]';
Pk = 1000 * Q;

Pest_kf = zeros(Ntest,2);
Pest_kf(1,:) = Xk(1:2);

for np = 2 : Ntest
   
    % Prediction
    X_k1_k = F * Xk;
    P_k1_k = F * Pk * F' + Q;
    
    % Measurements Prediction:
    Zk = H * X_k1_k;
    Sk = H * P_k1_k * H' + R;
    
    % Update:
    G = (P_k1_k * H') * pinv(Sk);
    E_k1_k = Pest_lls(np,:) - Zk;
    X_k1_k1 = X_k1_k + G * E_k1_k;
    P_k1_k1 = (eye(Xs) - G * H) * P_k1_k;
    
    Pest_kf(np,:) = X_k1_k1(1:2);
    
    Xk = X_k1_k1;
    Pk = P_k1_k1;
    
end

figure; grid on; hold on
plot(Pest_kf(:,1), Pest_kf(:,2), '-o');
plot(Pest_lls(:,1), Pest_lls(:,2), 'x');

%% Extended Kalman Filter:
close all
% CA model in polar coordinates
f_CA = @(x,dt)[x(1) + x(4) * dt .* cos(x(3)) + .5 .*x(5) * (dt^2) .* cos(x(3)); ... X
               x(2) + x(4) * dt .* sin(x(3)) + .5 .*x(5) * (dt^2) .* sin(x(3)); ... Y
               x(3);                ... HEADING
               x(4) + x(5) * dt;    ... SPEED
               x(5)];               ... ACCELERATION
% Derivative of CA model (Jacobian)
df_CA=@(x,dt)[1, 0, (-x(4) .* sin(x(3)) * dt - .5 * (dt^2) * x(5) .* sin(x(3))), (cos(x(3)) * dt),(.5 * (dt^2) * cos(x(3)));...
              0, 1, (x(4) .* cos(x(3)) * dt + .5 * x(5) * (dt^2) .* cos(x(3))), (sin(x(3)) * dt),(.5 * (dt^2) * sin(x(3)));...
              0, 0, 1, 0, 0; ...
              0, 0, 0, 1, dt;...
              0, 0, 0, 0, 1];
          
T = 1;             % measurments update 1 s;
Xs = 5;
Q = 1e-4 * diag([4e2 1e0 4e0 1e-9 1e-6]);

% Measurements model:
H = [1 0 0 0 0; ...
     0 1 0 0 0];
 
R = 1*diag([2 2]);

% Initialization:
Xk = [Pest_lls(1,:) pi/2 .1 .002]';
Pk = 20 * Q;

Pest_kf = zeros(Ntest,2);
Pest_kf(1,:) = Xk(1:2);

for np = 2 : Ntest
   
    % Prediction
    X_k1_k = f_CA(Xk, T); % Mean prediction
    
    F_CA = df_CA(X_k1_k, T);        % Jacobian derivation
    P_k1_k = F_CA * Pk * F_CA' + Q; % Covariance prediction
    
    % Measurements Prediction: (H is linear, no need to derive the Jacobian)
    Zk = H * X_k1_k;                % Mean 
    Sk = H * P_k1_k * H' + R;       % Covariance
    
    % Update:
    G = (P_k1_k * H') * pinv(Sk);
    E_k1_k = Pest_lls(np,:) - Zk;
    X_k1_k1 = X_k1_k + G * E_k1_k;
    P_k1_k1 = (eye(Xs) - G * H) * P_k1_k;
    
    Pest_kf(np,:) = X_k1_k1(1:2);
    
    Xk = X_k1_k1;
    Pk = P_k1_k1;
    
end

figure; grid on; hold on
plot(Pest_kf(:,1), Pest_kf(:,2), '-o');
plot(Pest_lls(:,1), Pest_lls(:,2), 'x');