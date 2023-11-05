%extra

% a = [ones(1, K); zeros(N*2, K)];
% for k = 1:K
%    x_train_batch = x_train(1:batch_size,:);
%    y_train_batch = y_train(1:batch_size, :);
%    for t = 1:batch_size       
%        for k1 = 1:K
%        
%            z(k1) = a(:,k1)'*x_train_batch(t,:);
% 
%            y_e(k) = z(k)/sum(z);
%        end
%        loss = loss + (y_train_batch - y_e(k))*x_train_batch(t,:);
%    end
%    a(:,k) = a(:,k) + mu*loss;
% end
% labels = zeros(K,1);
% for i = 1 : K
%     labels(i) = Class(track_s1(i,:));
% end
% 
% 
% z(k) = zeros(400,1);
% 
% z(k) = a(k)'*x + b(k);
% 
% for i = 1: 400
%     y_e(i) = sigma(a(i,:)*X_test(1, :)');
% end
% Pr = Pr./sum(Pr);
% y_e = 
% 
% %one got encoding
% lables_ohe = zeros(5000,400);
% %Y2 = zeros(5000,100);
% for i = 1:100
%     for j = 1+(i-1)*50:i*50
%         lables_ohe(j,labels(i)) = 1;
% 
%     end
% end






% OHE = []
% 
% 
% 
% 
% mu = 0.01;
% lamda = 5;
% N_e = 500;
% 
% 
% a = zeros(400,3);
% for e = 1:N_e%epochs
%     for l = 1:size(X_train, 1)
%         for i = 1 :size(a, 2)
%             a(Y_train(l), i) = a(Y_train(l), i) + mu*(Y_train(l) - sigma(a(Y_train(l),:)*X_train(l,:)'))*X_train(l,i); %;- mu*lamda*a(1, i) + mu*(Y_train(l) - sigma(a(:,1)'*X_train(l,:)'))*X_train(l,1);
%         end
%     end
% end
% 
% Pr = zeros(400,1);
% for i = 1: 400
%     Pr(i) = sigma(a(i,:)*X_test(1, :)');
% end
% Pr = Pr./sum(Pr);
% 
% %res = softmax(a, X_test(1, :))