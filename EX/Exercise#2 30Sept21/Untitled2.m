% load('ADSP_exercising30Sept21.mat')
% plot(t,y);
% z = y(4430:6800);
% roots = roots(z)
%zplzne(roots, 1)

y_trunc=y(4430:6800);
h=[1, -1];
xp=conv(y_trunc,h);
sigma0=1;
w=sigma0*randn(size(xp));
x=xp+w;

H = convmtx(y_trunc, 2);
h_e0 = pinv(H)*x;


x_e0 = conv(h_e0, y_trunc);
eps0 = (x - x_e0);
MSE0 = eps0.'* eps0/length(eps0);

i = 0;
X = zeros(2372, 20);
sigma = (0:0.5:10);

for i= 1:21
    %for i = 1:2372
        X = xp(i) + sigma(i) * randn(size(xp));
        h_e = pinv(H)*X;
        x_e = conv(h_e, y_trunc);
        eps = (x - x_e);
        MSE(i) = eps.'* eps/length(eps);
        
   % end

end

plot(sigma,log(MSE))
grid on
xlabel('sigma') 
ylabel('MSE, dB') 