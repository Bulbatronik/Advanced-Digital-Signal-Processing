%% Vibration data from Hossein Safarzadeh
load DataVibrations
data=table2array(sample);
channel=3;
x=data(:,channel);

%% Generation sinusoid+noise
Nt=10^6; 
t=[1:Nt]';
omega_o=2*pi*6/17;
sigma_w=.1; % max 50 for visualization, Nfft=2^15 or 2^16
x=cos(8*pi/17*t+2*pi*rand)+sigma_w*randn(Nt,1);

%% Periodogram
Nfft=1024;
% window=boxcar(Nfft);  % rect window
window=bartlett(Nfft);  % triangular window
Noverlap=Nfft/2;
% [Pxx,omega] = periodogram(x); % plain periodogram
[Pxx,omega] = pwelch(x,window,Noverlap,Nfft);
    figure(1)
    plot(omega,10*log10(Pxx),'k')
    % plot(omega,Pxx,'k')
    % plot(omega,pi*Pxx,'--k',omega,sigma_w^2*ones(size(omega)),'r'); % pi scaling to compensate for PSD in rad
    xlabel('\omega'); ylabel('PSD_x [dB]');

%% AR Spectral Analysis
N_AR=20;
[Pxx,omega] =  pyulear(x,N_AR,512,2*pi);
    figure(1)
    plot(omega,10+10*log10(Pxx),'-r')
    % plot(omega,pi*Pxx,'k') % ,omega,sigma_w^2*ones(size(omega)),'r'); % pi scaling to compensate for PSD in rad
    xlabel('\omega'); ylabel('PSD_x [dB]')

%% AR Pole plot
AR_mod=aryule(x,N_AR);
    figure(2)
    zplane([1],AR_mod);
    
%% Suggestions
% - check periodogram with datavibrations and then find the N_AR for AR to have same spectral-peaks (resonances)
% - check periodogram for sinusoid with large sigma_w (say sigma_w>=50 and
% find the length of data segmentation (Nfft) till there is a visibility of
% the sinusoid line
% - verify the distribution of poles for AR of sinusoid in noise for
% varying sigma_w
% - modify the code and test different options and other spectral analysis
% tools (search for spectral analysis in help line)
%                                   enjoy :) US 05Dec2018

