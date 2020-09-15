



close all;
clear all;

std_ = [2;2];
std2_Mat = diag(std_).^2;
% std2_Mat = [32,15;15,40];
% std_ = sqrt(std2_Mat);
l = (std_).*randn(2,1e4);
mean1 = mean(l');
figure;plot(l(1,:),l(2,:),'.');
hold on;
plot(mean1(1),mean1(2),'ro')

%%
f = @(x)([x(1)+x(2);0.1*x(1)^2+x(2)^2]);
% f = @(x)([x(1)+x(2);0.1 * x(1)^2 + 10 * x(2)^3]);
% f = @(x)([x(1)+x(2); 2* x(1) + 20* x(2)]);

for i = 1:length(l)
    l2(:,i) = f(l(:,i));
end
figure;plot(l2(1,:),l2(2,:),'.');
axis equal
mean(l2')
std_2 = std(l2');

% f(mean1)
sigma0 = mean1;

%% ukf mean
dim = 2;
lamada = 2;         % lamada的调整在很大程度上并不会改变对均值估计的结果
coeff = sqrt(lamada + dim);
sigma_ = getsigma(sigma0 , coeff*std_);
w0m = lamada/(dim+lamada);
w = 1/(dim*(dim+lamada));
sum_ = w0m*f(sigma0);
for i = 1:length(sigma_)
    sum_ = sum_ + w * f(sigma_(i,:));
end
mean_prior = sum_

list = l2-mean_prior;
cov_ = sqrt(list(1,:)*list(2,:)'/length(list));
[std_2(1) ,cov_;
    cov_ , std_2(1)]
%% 方差
alpha = 1;
beta = 2;
k = 3-2;
w0c = lamada/(dim + lamada) + 1 - alpha^2 + beta;
% w0c = 0.4;

sum_ = w0c*( (f(sigma0)- mean_prior)*(f(sigma0)- mean_prior)');
for i = 1:length(sigma_)
    sigmatp = sigma_(i,:);
    sum_ = sum_ + w * ( (f(sigmatp)- mean_prior)*(f(sigmatp)- mean_prior)');
end
std2_prior = sum_;
% std2_prior(4) = std2_prior(4)/2.7;
sqrt(std2_prior)


l3 = mean_prior + sqrt(std2_prior)*randn(dim,1e4);
hold on;
plot(l3(1,:),l3(2,:),'.')
legend 1 2
%% 



function out = getsigma(mean_,std_)
out(1,:) =  mean_+ [std_(1) , 0];
out(2,:) =  mean_+ [-std_(1) , 0];
out(3,:) =  mean_+ [ 0 , std_(2)];
out(4,:) =  mean_+ [0, -std_(2)];
end

