clear all
N=128*1024;
omega=13*pi/64; %be aware on bias!
A=1;
M=40;

% SNR=20;
% sigma_w=A/sqrt(2*SNR);

sigma_w=100; 

S_est=zeros(N,1);
for m=1:M,
    % data generation
    x=A*cos(omega*[0:N-1]'+2*pi*rand)+sigma_w*randn(N,1);
    % SIMPLE periodogram    
    P=(abs(fft(x)).^2)/N;
    S_est(:,1)=S_est(:,1)+P;
end
S_est=S_est/M;

semilogy(2*pi*[-N/2:N/2-1]/N,fftshift(S_est),'-k',omega,10,'+')
xlabel('\omega')
ylabel('S_e_s_t(\omega)')
grid on
    

