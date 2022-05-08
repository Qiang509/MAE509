%% HW1 Question 6-a

%% MATLAB script for
% scatterplots of the project
%     gamma = sqrt((2*(a+b)*h)/(a*b*k))
%     c1 = -(Q/(k*gamma))*((exp(gamma*L)*(h+k*gamma))/(exp(-gamma*L)*(h-k*gamma)+exp(gamma*L)*(h+k*gamma))
%     c2 = Q/(k*gamma) + c1
%     T = c1*exp(-gamma*x)+c2*exp(gamma*x)+T_amb
	
clear all; close all; clc	

%  coefficients	

a = 0.95;      % cm
b = 0.95;      % cm
L = 70.0;      % cm
k = 2.37;      % W/cm C
T_amb = 21.29; % C

% number of samples
N = 5000; %tested 1k,10k,100k

% drawing samples from parameters and construct sample matrix M
theta1 = unifrnd(-36,0,[N,1]); % Q
theta2 = unifrnd(.001,.003,[N,1]);	% h
x = .1;
M = [theta1, theta2];

% compute vector of model outputs
gamm = sqrt((2*(a+b).*theta2)./(a*b*k));
c1 = -(theta1./(k*gamm)).*((exp(gamm*L).*(theta2+k.*gamm))./(exp(-gamm.*L).*(theta2-k.*gamm)+exp(gamm*L).*(theta2+k.*gamm)));
c2 = theta1./(k*gamm) + c1;
T = c1.*exp(-gamm.*x)+c2.*exp(gamm.*x)+T_amb;

%% scatter plots
figure
plot(M(:,1), T, '*b')
xlim([-36 0]), grid on
axis square,xlabel('\theta_1'),ylabel('y')		
set(gca,'FontSize',24)	
ylim([20 400])
print('x1y','-dpng')	
	
figure
plot(M(:,2), T, '*b')
xlim([.001 .003]), grid on		
axis square,xlabel('\theta_2'),ylabel('y')	
ylim([20 400])
set(gca,'FontSize',24)	
print('x2y','-dpng')

figure
plot(M(:,1),M(:,2),'*')
