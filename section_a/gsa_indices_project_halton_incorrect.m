%% Project

%%
%  This code illustrates the implementation of the 
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

clear all, clc
close all

%% Setup the model and define input ranges
%  coefficients

x  = [50];

% number of parameters
p = 2;

% parameter ranges
param1 =  [-36 0];
param2 =  [.001 .003];

%% Sample parameter space:
% number of samples
M = 100000

% Compute [A], [B] and [C] as random variables
A(:,1) = param1(1) + (param1(2) - param1(1)).*net(haltonset(1,'Skip',100),M);
A(:,2) = param2(1) + (param2(2) - param2(1)).*net(haltonset(1,'Leap',3),M);

B(:,1) = param1(1) + (param1(2) - param1(1)).*net(haltonset(1,'Skip',30),M);
B(:,2) = param2(1) + (param2(2) - param2(1)).*net(haltonset(1),M);

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
VARy = mean(yA.^2) - f0^2;

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
(STi(1)-Si(1))+(Si(1)+Si(2))

img1 = figure(1)
bar(abs(indices))
xlabel('\theta'),ylabel('Sensitivity Indices'), grid on		
set(gca,'FontSize',24)
legend('first-order', 'total effects','Location','bestoutside')
title(['x=' num2str(x)])

% scatter plots
img2 = figure(2)
plot(A(:,1), yA, '*b')
xlabel('\phi [W/m^{2}]'),ylabel('Y [C°]'), grid on	
ylim([0 160]), xlim(param1)
set(gca,'FontSize',24)
title(['x=' num2str(x)])
	
img3 = figure(3)
plot(A(:,2), yA, '*b')
xlabel('h [W/(m^{2}*K)]'),ylabel('Y [C°]'), grid on
ylim([0 160]), xlim(param2)
set(gca,'FontSize',24)	
title(['x=' num2str(x)])	

figure(4)
plot(A(:,1),A(:,2),'*')
xlabel('\theta_1'), ylabel('\theta_2'), grid on
set(gca,'FontSize',24)
title('h vs \phi')