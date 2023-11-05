clear all; close all; clc

%% Frequency Estimation, Ch. 9

A = 1;  % amplitude
N = 128; % # of time samples
Oversampling = 4; 
M = N * Oversampling;  % # of FFT points
n = ( 0 : N - 1 )';   % time axis
SNR = - 30 : 45;
snr = 10 .^ ( SNR / 10 );  % linear snr
Pn = A ^ 2 / 2 ./ snr;  % power of the noise
Nrun = 1e3; % # of montecarlo sims
fo = 0.2;  % normalized frequency

% bounds
crb = 12 / N / ( N ^ 2 - 1 ) / 4 / pi ^ 2 ./ snr;  % crb for freq estimation
MSE_floor = 1 / 12 / M ^ 2;  % floor for high snr
MSE_ceil = 1 / 4 / 12;  % ceil for low snr

% init variables
f_est1 = zeros( Nrun , 1 ); 
f_est2 = zeros( Nrun , 1 );

MSE1 = zeros( size( snr ) );
MSE2 = zeros( size( snr ) );

for isnr = 1 : length( snr )
    
    for irun = 1 : Nrun
        
        phi = rand * 2 * pi;  % phi~U(0,2*pi)
        w = sqrt( Pn( isnr ) ) * randn( N , 1 );  % w ~ N(0,Pn)
        x = A * cos( 2 * pi * fo .* n + phi ) + w;  % rx signal
        
        %estimate freq
        X = fft( x , M );  % M-points fft
        S = abs( X ) .^ 2 / N; % periodogram
        
        [ ~ , f_est1( irun ) ] = max( S( 2 : M / 2 ) ); % max of periodogram, positive freq (zero exluded)
        
        % parabolic interpolation
        f_cent = f_est1( irun ) + 1; % center
        
        num = S( f_cent - 1 ) - S( f_cent + 1 );
        den = 2 * ( S( f_cent - 1 ) + S( f_cent + 1 ) - 2 * S( f_cent ) );
        
        df = num / den; % correction
        
        f_est2( irun ) = f_est1( irun ) + df; % freq estimate by parabolic interp
        
    end
    
    MSE1( isnr ) = mean( ( f_est1 / M - fo ) .^ 2 );  
    MSE2( isnr ) = mean( ( f_est2 / M - fo ) .^ 2 );
    
end

figure( 1 )
semilogy( SNR , MSE1 , '-' , SNR , MSE2 , '-o' , SNR , crb , SNR , MSE_floor * ones( size( snr ) ) , ':' ,...
    SNR , MSE_ceil * ones( size( snr ) ) , ':' ,'linewidth' , 2 , 'markersize' , 8 );
leg = legend( 'MSE' , 'MSE_{parbolic}' , 'CRB' , 'MSE_{floor}' , 'MSE_{ceiling}' );
set( leg , 'fontsize' , 14 )
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
        
        
        
        
        
        
        
        
