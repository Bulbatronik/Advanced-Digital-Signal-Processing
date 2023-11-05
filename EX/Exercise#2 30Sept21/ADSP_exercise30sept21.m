% Exercising ADSP on 30Sept2021
% y from ADSP_exercising30Sept21.mat is voice (but can be anything)
% x'=y*h and x=y*h+w

clear all % clear all variables before starting
%%
load ADSP_exercising30Sept21
% plot of roots of segmented voice file over meaningful section
y=y(4430:6800); %segmentation and reassignment to y
sound(y,8E3)
r=roots(y);

plot(r, 'o'); axis('square'); 
%%
sigma_w=1;

h=[1, -1];
x_con=conv(y,h);

w=sigma_w*randn(size(x_con));
x=x_con+w; % optional 
Y=zeros(length(y)+length(h)-1,length(h));
for i=1: length(h)
Y(i:length(y)-1+i,i)=y;
end
h_est=pinv(Y)*x;
%% loop to compute the MSE vs sigma
sigma_w=logspace(-1,2,100); %set of noise power

MSE_h=zeros(size(sigma_w));
MSE_x=zeros(size(sigma_w));
MSE_xclean=zeros(size(sigma_w));

h=[1, -1]';
x_con=conv(y,h);

for k_sigma=1:length(sigma_w),
    w=sigma_w(k_sigma)*randn(size(x_con));
    x=x_con+w; % optional
    Y=zeros(length(y)+length(h)-1,length(h));
    for i=1: length(h)
        Y(i:length(y)-1+i,i)=y; % convolution matrix Y
    end
    h_est=pinv(Y)*x;
    MSE_h(k_sigma)=mean((h_est-h).^2); %MSE of filters' tap
    MSE_x(k_sigma)=mean((Y*h_est-x).^2)/length(x); %MSE of data-out (practical choice in many context)
    MSE_xclean(k_sigma)=mean((Y*h_est-x_con).^2)/length(x); % MSE of data out of filter (specific analysis on filer)
end

figure(1)
loglog(sigma_w.^2,MSE_h)
title('MSE_h vs \sigma^2_w')
xlabel('\sigma_w^2')
ylabel('MSE_h')

figure(2)
loglog(sigma_w.^2,MSE_x)
title('MSE_x vs \sigma^2_w')
xlabel('\sigma_w^2')
ylabel('MSE_x')

figure(3)
loglog(sigma_w.^2,MSE_xclean)
title('MSE_x_'' vs \sigma^2_w')
xlabel('\sigma_w^2')
ylabel('MSE_xclean')

%%
%pause

%% Alternative plots
% figure(4)
% subplot(221)
% zplane(y',[1]); % built-in function to plot poles and zeros of LTI systems
% title('poles and zeros of y')
% 
% subplot(222)
% loglog(sigma_w.^2,MSE_h)
% grid('on')
% title('MSE_h vs \sigma^2_w')
% xlabel('\sigma_w^2')
% ylabel('MSE_h')
% 
% subplot(223)
% loglog(sigma_w.^2,MSE_x)
% grid('on')
% title('MSE_x vs \sigma^2_w')
% xlabel('\sigma_w^2')
% ylabel('MSE_x')
% 
% subplot(224)
% loglog(sigma_w.^2,MSE_xclean)
% grid('on')
% title('MSE_x_'' vs \sigma^2_w')
% xlabel('\sigma_w^2')
% ylabel('MSE_xclean')
% 
% 
% 








