%% Audio recorder:

clear, clc, close all;

T = 1;      % Record duration in seconds
Fs = 8e3;   % sample rate 


recObj = audiorecorder;

tic
fprintf('Recording starts in 3 seconds [ ')
while toc < 2
    fprintf('. ')
    pause(0.99)
end
fprintf(']\n');
disp('Start speaking.')
pause(0.3)
recordblocking(recObj, T);
disp('End of Recording.');

% Uncomment to play the recorder audio:
% pause(2);
% play(recObj);

y = getaudiodata(recObj);
t = (0:length(y)-1)/Fs;
N = (1:length(t));
% plot the audio
figure;
plot(t, y)%
figure;
stem(N, y);
%%
Y=fft(y);
