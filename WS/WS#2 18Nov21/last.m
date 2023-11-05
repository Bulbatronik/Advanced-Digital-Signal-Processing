clear;
clc
close all
load data_1.mat

v = 335;%m/s

phi=1:800;
r_r = sqrt(25-2*cosd(phi)+0.04);
r_l = sqrt(25+2*cosd(phi)+0.04);

t_r = sqrt(25-2*cosd(phi)+0.04)/v;
t_l = sqrt(25+2*cosd(phi)+0.04)/v;

time_s = 1/Fs * length(s);%duration of the audio

num_sample = 360;% # of samples
time_g = time_s/num_sample;%duration of each sample

sample_size = round(length(s)/num_sample); %size of each sample

x_l_sampled = zeros(sample_size,num_sample);
x_r_sampled = zeros(sample_size,num_sample);

l_r_cross = zeros(sample_size*2-1,num_sample);%2 times the signal

k=1;
%each column represents one sample
for i=0:num_sample-1
    ind = i*sample_size+1;
    if(i==num_sample-1)
        x_l_sampled(1:length(s)-ind+1,i+1)=s(ind:end,1);
        x_r_sampled(1:length(s)-ind+1,i+1)=s(ind:end,2);
    else
        x_l_sampled(:,i+1)=s(ind:(i+1)*sample_size,1);
        x_r_sampled(:,i+1)=s(ind:(i+1)*sample_size,2);
    end
    l_r_cross(:,i+1) = xcorr(x_l_sampled(:,i+1),x_r_sampled(:,i+1));
    [values(k),indexes(k)] = max(l_r_cross(:,i+1));
    indexes(k) = indexes(k) - sample_size;
    k=k+1;
end

indexes = indexes/Fs;

lags = -sample_size+1:sample_size-1;

for i=1:num_sample

    plot(lags,l_r_cross(:,i))
end
title('Left-Right cross correlation')

%%
%graph
%question 1
figure
subplot(2,1,1)
plot(r_l)
xlabel('degree')
ylabel('distance (m)')
title('left distance')
subplot(2,1,2)
plot(r_r)
title('right distance')
xlabel('degree')
ylabel('distance (m)')

figure
subplot(2,1,1)
plot(t_l)
xlabel('degree')
ylabel('time (s)')
title('left time')
subplot(2,1,2)
plot(t_r)
title('right time')
xlabel('degree')
ylabel('time (s)')

figure
plot(t_r-t_l)
title('time difference (s)')
ylabel('difference')
xlabel('degree')

plot(indexes)