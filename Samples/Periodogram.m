% Matlab Code used for DSP Course class exercise 08Jan.2014
% Similar problems can be topic for oral exams
% U.Spagnolini
%
clear all

% filter design
[b,a]=cheby1(15,.1,.4/2);
w=pi*linspace(-1,1,201);
h=freqz(b,a,w);
% semilogy(w,abs(h))

% simulation parameters
Ntot=1000;
N=500;
M=Ntot/N;
f=[-(N/2):N/2-1]/(N/2); % frequency span for plot purposes

n=randn(Ntot,1);
x=filter(b,a,n);
x=x+.05*cos(pi*.25*[1:Ntot]');

win=bartlett(N); % boxcar=no-window!

for m=1:M,
    X(m,:)=fft(win.*x(1+N*(m-1):N+N*(m-1)),N)';
    X(m,:)=fftshift(X(m,:));
end
% plot(f,(1/N)*abs(X).^2,'-',w/pi,abs(h).^2,'--')
semilogy(f,(1/N)*mean(abs(X).^2),'-',w/pi,abs(h).^2,'--')
axis([-1,1,10^-6,10])
xlabel('frequency (normalized)')
ylabel('PDS')
title('Solid: WOSA, dashed: true PSD')


