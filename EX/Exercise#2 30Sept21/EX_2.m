% load('ADSP_exercising30Sept21.mat')
% plot(t,y);
% z = y(4430:6800);
% roots = roots(z)
%zplzne(roots, 1)
clear all
load ADSP_exercising30Sept21
y_trunc=y(4430:6800);
h=[1, -1];
xp=conv(y_trunc,h);
sigma0=1;
w=sigma0*randn(size(xp));
x=xp+w;

Y = convmtx(y_trunc, 2);
h_e0 = pinv(Y)*x;

% for i=1: length(h)
%     Y(i:length(y)-1+i,i)=y; % convolution matrix Y
% end

x_e0 = conv(h_e0, y_trunc);
eps0 = (x - x_e0);
MSE0 = eps0.'* eps0/length(eps0);

i = 0;
X = zeros(2372, 20);
sigma =logspace(-1,2,100);

for i= 1:length(sigma)
    %for i = 1:2372
        X = xp(i) + sigma(i) * randn(size(xp));
        h_e = pinv(Y)*X;
        x_e = conv(h_e, y_trunc);
        eps_x = (X - x_e);
        eps_h = (h_e - h.');
        %MSE(i)=(x_e - X).^2/length(x);
        MSE_x(i) = eps_x.'* eps_x/length(eps_x);
        MSE_h(i) = eps_h.'* eps_h/length(eps_h); %MSE of filters' tap
        %MSE_xclean(k_sigma)=mean((Y*h_est-x_con).^2)/length(x); % MSE of data out of filter (s
   % end

end

figure(1)
loglog(sigma.^2,MSE_h)
xlabel('sigma') 
ylabel('MSE_h, dB') 

figure(2)
loglog(sigma.^2,MSE_x)
xlabel('sigma') 
ylabel('MSE_x, dB') 