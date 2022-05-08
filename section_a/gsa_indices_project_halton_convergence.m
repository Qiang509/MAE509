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
VARyv = []; Mv = []; S1v = []; S1Tv = []; S2v = []; S2Tv = [];
x  = [50];

% number of parameters
p = 2;

% parameter ranges
param1 =  [-36 0];
param2 =  [.001 .003];

%% Sample parameter space:
% number of samples
for M = [1:2:100]*10^3;
A = []; B = []; C = [];
Mv = [Mv, M];

halt = net(haltonset(4),M);

% Compute [A], [B] and [C] as random variables
A(:,1) = param1(1) + (param1(2) - param1(1)).*halt(:,1);
A(:,2) = param2(1) + (param2(2) - param2(1)).*halt(:,2);

B(:,1) = param1(1) + (param1(2) - param1(1)).*halt(:,3);
B(:,2) = param2(1) + (param2(2) - param2(1)).*halt(:,4);


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
VARyv = [VARyv VARy];

for i = 1:p
    yCi = yC(:,i);

	% first order indices	
    Si(i)  = ( 1/M*sum(yA.*yCi) - f0^2 ) / VARy ; 

    % total effects indices
    STi(i) = 1 -  ( 1/M*sum(yB.*yCi) - f0^2 ) / VARy ;
end

%% Plot results
% sensitivity indices
indices = [Si' STi'];
S1v = [S1v Si(1)]; S1Tv = [S1Tv STi(1)]; S2v = [S2v Si(2)]; S2Tv = [S2Tv STi(2)];

% img1 = figure(1)
% bar(abs(indices))
% ylim([0 1])
% xlabel('\theta'),ylabel('T = c1*exp(-gamma*x)+c2*exp(gamma*x)+T_amb'), grid on		
% set(gca,'FontSize',24)
% legend('first-order', 'total effects')
% title(['x=' num2str(x)])
% 
% % scatter plots
% img2 = figure(2)
% plot(A(:,1), yA, '*b')
% xlabel('\theta_1'),ylabel('Y'), grid on	
% ylim([0 160]), xlim(param1)
% set(gca,'FontSize',24)
% title(['x=' num2str(x)])
% 	
% img3 = figure(3)
% plot(A(:,2), yA, '*b')
% xlabel('\theta_2'),ylabel('Y'), grid on
% ylim([0 160]), xlim(param2)
% set(gca,'FontSize',24)	
% title(['x=' num2str(x)])	
    if M == 9000
        figure(3)
        plot(A(:,1),A(:,2),'*')
        title('Halton \phi vs h for M = 9000')
        xlabel('\phi'), ylabel('h')
    end
end

figure(1)
plot(Mv,S1v)
hold on 
plot(Mv,S2v)
xlabel('# of Samples')
ylabel('Total Indices')
legend('S_{1}','S_{2}')
title('Halton Convergence of Indices')


figure(2)
plot(Mv,VARyv)
xlabel('# of Samples')
ylabel('Variance of Output')
title('Halton Convergence of Variance in Output')

% disp(VARyv)