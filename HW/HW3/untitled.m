%HW3 Google
clc, clear, close all;

load('GoogleDataset.mat');

K = 20;
T20 = 397;
%% 

full_click= xcorr(Click(1:396, 1));
Corr_click = full_click(396:end);

%x_est = zeros(396,396);
%i=1;
%R = toeplitz(Correlations_click);


for step = 1:396
    for length = 1:396
        R = toeplitz(Corr_click(1:length));
        p =  - Corr_click(step:step+length);




        a = -inv(R)*p;
    


   
%         x_est = zeros(396);
%         while(n<397)
%             x_est(n) = a'* Click(n-length(a),1);     
%             n=n+1;
     end
end

figure;
    plot(Conversion(1:396,1))
    hold on
    plot(x_est)


% heatmap(abs(Correlations_click));
% 
% 
% 
% ind = 1:19;
% coefs =  ([ind; sum(abs(Correlations_click))])';
% ranked = sortrows(coefs, 2);


%% 
% figure;plot(Conversion(1:396,1:19))
% title('Conversion')
% xlabel('time')
% ylabel('Conversion')
% lgd = legend(string(1:19));
% title(lgd,'company')
% 
% 
% figure;plot(Cost(1:396,1:19))
% title('Cost')
% xlabel('time')
% ylabel('Cost')
% lgd = legend(string(1:19));
% title(lgd,'company')
% 
% figure;plot(Click(1:396,1:19))
% title('Click')
% xlabel('time')
% ylabel('Click')
% lgd = legend(string(1:19));
% title(lgd,'company')





% 
% %Click prediction
% Correlations_click = corrcoef(Click(1:397, 1:19));
% heatmap(abs(Correlations_click));
% 
% ind = 1:19;
% coefs =  ([ind; sum(abs(Correlations_click))])';
% ranked = sortrows(coefs, 2);
