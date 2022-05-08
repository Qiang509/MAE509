%% Project

%% Saltelli estimators of Sobol indices
%  This code illustrates the implementation of the Monte Carlo estimators
%  for computing the first first-order indices 
%  and total effects indices 
%
%     T = c1*exp(-gamma*x)+c2*exp(gamma*x)+T_amb
%     gamma = sqrt((2*(a+b)*h)/(a*b*k))
%     c1 = -(Q/(k*gamma))*((exp(gamma*L)*(h+k*gamma))/(exp(-gamma*L)*(h-k*gamma)+exp(gamma*L)*(h+k*gamma))
%     c2 = Q/(k*gamma) + c1
%
%     parameters theta = [theta1, theta2]
%     theta1 = Q; theta2 = h
%     x held constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc
close all;

%% Setup the model and define input ranges
%  coefficients
for x  = [10:4:66]
% number of parameters
p = 2;
% parameter ranges
param1 =  [-21 -15];
param2 =  [.00191-(3e-4) .00191+(3e-4)];
%param3 = [10 66];

%% Sample parameter space:
% number of samples
M = 10000; %tested 25k, 50k, 100k 

% Compute [A] and [B] as random variables
A(:,1) = param1(1) + (param1(2) - param1(1)).*lhsdesign(M,1);
A(:,2) = param2(1) + (param2(2) - param2(1)).*lhsdesign(M,1);
%A(:,3) = param3(1) + (param3(2) - param3(1)).*linspace(0,1,M);

B(:,1) = param1(1) + (param1(2) - param1(1)).*lhsdesign(M,1);
B(:,2) = param2(1) + (param2(2) - param2(1)).*lhsdesign(M,1);
%B(:,3) = param3(1) + (param3(2) - param3(1)).*linspace(0,1,M);

C = zeros(M,p,p);
for i = 1:p
    C(:,:,i) = B;
    C(:,i,i) = A(:,i);
end

%% Run the model and compute selected model output at sampled parameter
for  j = 1:M
    yA(j,1) = project_ind(A(j,:),x);
    yB(j,1) = project_ind(B(j,:),x);
    for i = 1:p
        yC(j,i) = project_ind(C(j,:,i),x);
    end
end

%% Compute sensitivity indices
f0  = mean(yA) ;
VARy = mean(yA.^2) - f0^2 ;

for i = 1:p
    yCi = yC(:,i);

	% fist order indices	
    Si(i)  = ( 1/M*sum(yA.*yCi) - f0^2 ) / VARy ; 

    % total effects indices
    STi(i) = 1 -  ( 1/M*sum(yB.*yCi) - f0^2 ) / VARy ;
end

%% Plot results
% sensitivity indices
indices = [Si' STi']

figure
bar(abs(indices))
axis square,xlabel('\theta'),ylabel('Y = sin(\theta_1) + a*sin(\theta_2)^2 + b*\theta_3^4*sin(\theta_1)'), grid on		
set(gca,'FontSize',24)
legend('first-order', 'total effects')

% scatter plots
figure
plot(A(:,1), yA, '*b')
axis square,xlabel('\theta_1'),ylabel('Y'), grid on		
set(gca,'FontSize',24)		
	
figure
plot(A(:,2), yA, '*b')
axis square,xlabel('\theta_2'),ylabel('Y'), grid on		
set(gca,'FontSize',24)	

% figure
% plot(A(:,3), yA, '*b')
% axis square,xlabel('x'),ylabel('Y'), grid on		
% set(gca,'FontSize',24)	
end