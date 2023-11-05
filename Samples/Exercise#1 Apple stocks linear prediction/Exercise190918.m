n=500;
m=600;
p=20;
clear y
clear X

for k=n:n+m,
    xrow=[1,x(k-p+1:k-1)'];
    X(k-n+1,:)=xrow;
    y(k-n+1)=x(k);
end
 y=y';

 theta=inv(X'*X)*X'*y;
 
 P=X*inv(X'*X)*X';
 
 rmse=sqrt((y'*y-y'*P*y)/m)
 
 %comparison between real and predicted
 figure
 plot(t(n:n+m),y,'-k',t(n:n+m),X*theta,'-r')
 xlabel('date-day'); title('values (black-line) vs predicted (red-line)')
 
 % histogram of absolute value of the error
 figure
 error_vector=y-P*y;
 hist(error_vector,20)
 
 % histogram of % value of the error (investement-oriented metric)
 figure
 error_vector_percentage=(y-P*y)./y;
 hist(100*error_vector_percentage,20)
 xlabel('% or error')
 
 % Matrix property: eigenvalues of projection matrix P are =0,1
figure
[Q,L]=eig(P);
plot(abs(diag(L)))