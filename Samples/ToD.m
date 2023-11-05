clear all;
close all;
clc

%% Time of Delay Estimation

N = 256; % # of samples per interval
M = 5e3; % intervals

Ntot = N * M;  % # tot number of samples

SNR = - 10 : 5 : 40;  % snr in dB
snr = 10 .^ ( SNR / 10 );  % snr linear

tau_vec = 0 : .1 : 1;  % values of tau

% signal generation
w = randn( Ntot , 1 );   % random noise
filter_order = 8;        % filter order
Bg = 0.05;               % positive bandwidth
[ b , a ] = butter( filter_order , 2 * Bg );   % butterworth filter
g = filter( b , a , w );   % noise filtered

% no delay
g1 = g;
G = fft( g );  % FFT of g
freq = [ 0 : Ntot / 2 , - Ntot / 2 + 1 : - 1 ]';  % freq = fftshift( [ - Ntot / 2 + 1 : Ntot / 2 ] );

% initi mse
mse = zeros( length( snr ) , length( tau_vec ) );

for itau = 1 : length( tau_vec )
    
    tau = tau_vec( itau ); % select value of tau
    
    % add delay in frequency domain
    G2 = G .* exp( 1i * 2 * pi * tau .* freq / Ntot );
    
    % time domain signal
    g2 = real( ifft( G2 ) );
    
    % noise for the overall signal duration
    w1 = randn( Ntot , 1 );
    w2 = randn( Ntot , 1 );
    
    for isnr = 1 : length( snr )
        
        % signals g1 and g2 + noise
        x1 = sqrt( snr( isnr ) ) * g1 + w1;
        x2 = sqrt( snr( isnr ) ) * g2 + w2;
        
        tau_est = zeros( M , 1 );
        
        for iM = 1 : M
            
            % indexes of iM block
            idx_m = ( iM - 1 ) * N + 1 : iM * N;
            
            % iM block of x1 and x2
            x1m = x1( idx_m );
            x2m = x2( idx_m );
            
            % find the crosscorr
            c = xcorr( x1m , x2m );
            
            % take the maximum
            [ ~ , idx_max ] = max( c );
            
            % parabolic interpolation
            num = c( idx_max - 1 ) - c( idx_max + 1 );
            den = 2 * ( c( idx_max - 1 ) + c( idx_max + 1 )...
                - 2 * c( idx_max ) );
            
            % correction
            p = num / den;
            
            % estimate tau
            tau_est( iM ) = idx_max - N + p;
            
        end
        % compute the error
        err = tau_est - tau;
        
        % compute the mse
        mse( isnr , itau ) = mse( isnr , itau ) + err' * err / M;
        
        fprintf( 'Completed iteration SNR = %d, tau = %d \n' , SNR( isnr ) , tau ) 
        
    end
    
end

Bg_eff = 1.1 * Bg;
crb = 3 / ( 8 * N * pi ^ 2 * Bg_eff ^ 3 ) * ( 1 + 2 .* snr ) ./ snr .^ 2;

% plot result
figure( 1 )
semilogy( SNR , sqrt( crb ) , '--' , SNR ,sqrt( mse ) , '-' )
xlabel( 'SNR [dB]' )
ylabel( 'RMSE' )
legend( 'CRB' )
        
    
    
    
    
    
    
    
    
    
    



