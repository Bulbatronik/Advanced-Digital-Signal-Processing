 clc;
 clear all;
 
 % Data File for marks of students in two subject 
 % to determine whether a student is pass or fail
 
 x=xlsread('marks.xlsx');
 ytrain=x(:,end); % Target variable
 xtrain=zscore(x(:,1:end-1));% Normalized Predictors
 
 p=length(x);


 xtrain=[ones(length(xtrain),1) xtrain]; % one is added for calculation of biases.
 xtest=xtrain;
 ytest=ytrain;
 
 %compute cost and gradient
  iter=1000; % No. of iterations for weight updation
  
  theta=zeros(size(xtrain,2),1); % Initial weights
  
  alpha=0.1 % Learning parameter
  
  [J grad h th]=cost(theta,xtrain,ytrain,alpha,iter) % Cost funtion
 
   ypred=xtest*th; %target prediction
 
   % probability calculation
 [hp]=sigmoid(ypred); % Hypothesis Function
 ypred(hp>=0.5)=1;
 ypred(hp<0.5)=0;
 
% Decision Boundary

syms x1 x2
fnn=th(1)+th(2)*x1+th(3)*x2-0.5

 figure
 hold on
 scatter(xtest(ytest==1,2),xtest(ytest==1,3),'b+','linewidth',5.0)
 scatter(xtest(ytest==0,2),xtest(ytest==0,3),'r','linewidth',5.0)
 h1=ezplot(fnn)
 set(h1,'color','r')
 legend('Pos class','Neg. class','Decision boundary')
 xlim([-2 2])
 ylim([-2 2])
 
