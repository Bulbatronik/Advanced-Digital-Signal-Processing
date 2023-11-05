%base predictor - average
%Conversion_average = Conversion(:,1:19)

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

%Click prediction
Correlations_click = corrcoef(Click(1:397, 1:19));
%heatmap(abs(Correlations_click));

ind = 1:19;
coefs =  ([ind; sum(abs(Correlations_click))])';
ranked = sortrows(coefs, 2);
%% 


%figure;plot(Click(:,1:19))