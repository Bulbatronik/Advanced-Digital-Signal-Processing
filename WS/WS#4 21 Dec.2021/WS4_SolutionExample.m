%% WS#4 - Array Processing - 21 Dec. 2021
clc, clear, close all;

%% 1. Trajectory Generation

dt = 0.5;

% c = [x, y, v, h], where x,y are the position, v is the speed and h is the
% direction of motion
f =@(c) [(c(1)+c(3)*cos(c(4))*dt), ...
         (c(2)+c(3)*sin(c(4))*dt), ... 
         c(3), ...
         c(4)];
     
Q = [0 0 0 0; ...
     0 0 0 0; ...
     0 0 0 0; ...
     0 0 0 pi/360] + abs(1e-3*diag(randn(4,1)));

K = 200;        % number of steps considered

% Generating a new trajectory
tr_1 = [40, 30, 2, 3*pi/2];     % initial state
for i = 2 : K
    
    tr_1(i,:) = f(tr_1(i-1,:)) + (sqrtm(Q) * randn(4,1)).';
    
end

figure; 
plot(tr_1(:,1), tr_1(:,2), '-x');
figure;
plot(tr_1(:,3))

%% 2. Optimal Beamformer

snr_db = 5:5:40;
N = 12;
v = 3*10^8;
fo = 2.9e9;
d = v / (2*fo);
lo = v / fo;
wo = 2*pi/lo;
pos_a = [zeros(N,1) (0:N-1).'*d];
L = 100;

% figure; hold on;
% plot(tr_1(:,1), tr_1(:,2), '-x');
% plot(pos_a(:,1), pos_a(:,2), 'o');

% compute the true distances and angles
rho_true = pdist2(pos_a(1,:), tr_1(:,1:2));
aoa_true = atan2(tr_1(:,2), tr_1(:,1));

% compute the true channel
h_true = exp(1i * wo * d * (0:N-1).' * sin(aoa_true.'));

az = -pi/2 : pi/180 : pi/2;
b = exp(1i * wo * d * (0:N-1).' * sin(az)) / sqrt(N);

sig2_s = 1;

for i = 1 : K
    
    h_c = h_true(:,i);
    
    for j = 1 : length(snr_db)
        snr_c = db2pow(snr_db(j));
        sig2_w = sig2_s / snr_c;
        
        w = sqrt(sig2_w/2) * (randn(L,1) + 1i* randn(L,1));
        s = sqrt(sig2_s/2) * (randn(L,1) + 1i* randn(L,1));
        
        y = h_c * s.' + w.';
        
        g(i,j,:) = mean(abs(b' * y).^2, 2);
        
    end
    
end

figure;
imagesc(squeeze(g(:,end,:)))
g_ideal = zeros(K,length(az));
for i = 1 : K
    g_ideal(i,abs(az-aoa_true(i)) == min(abs(az-aoa_true(i)))) = 1;
end
figure;
imagesc(g_ideal)




