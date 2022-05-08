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
v1 = VideoWriter('Halton_Sobol_Indices.mp4','MPEG-4');
v2 = VideoWriter('Halton_scattertheta_1.mp4','MPEG-4');
v3 = VideoWriter('Halton_scattertheta_2.mp4','MPEG-4');;
v1.FrameRate = 2; v2.FrameRate = 2; v3.FrameRate = 2;
open(v1); open(v2); open(v3); 
VARyv = []; xv = []; S1v = []; S1Tv = []; S2v = []; S2Tv = [];

for x  = [10:4:70];
xv = [xv, x];

% number of parameters
p = 2;
% parameter ranges
param1 =  [-36 0];
param2 =  [.001 .003];

%% Sample parameter space:
% number of samples
M =  100000;
halt = net(haltonset(4),M); % create 4 unique halton vectors

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
indices = [Si' STi']
S1v = [S1v Si(1)]; S1Tv = [S1Tv STi(1)]; S2v = [S2v Si(2)]; S2Tv = [S2Tv STi(2)];

img1 = figure(1)
bar(abs(indices))
ylim([0 1])
xlabel('\theta'), ylabel('Sensitivity'), grid on		
set(gca,'FontSize',24)
legend('first-order', 'total effects','Location','bestoutside')
title(['x=' num2str(x)])

% scatter plots
img2 = figure(2)
plot(A(:,1), yA, '*b')
xlabel('\phi [W/m^{2}]'),ylabel('Y [C°]'), grid on	
ylim([0 280]), xlim(param1)
set(gca,'FontSize',24)
title(['x=' num2str(x)])
	
img3 = figure(3)
plot(A(:,2), yA, '*b')
xlabel('h [W/(m^{2}*K)]'),ylabel('Y [C°]'), grid on
ylim([0 280]), xlim(param2)
set(gca,'FontSize',24)	
title(['x=' num2str(x)])	

% video writing
f1=getframe(img1);
writeVideo(v1,f1);
f2=getframe(img2);
writeVideo(v2,f2);
f3=getframe(img3);
writeVideo(v3,f3);
pause(0.1)

end

% theta1 vs theta 2
% figure(4)
% plot(A(:,1),A(:,2),'*')
% figure(5)
% plot(xv,VARyv)

% Sensitivity
figure(4)
plot(xv,S1v,'-o'), hold on
plot(xv,S2v,'-o')
plot(xv,(S1Tv-S1v)+(S1v+S2v),'-o') % e.q. (7)
plot(xv,(S1Tv-S1v),'-o')
ylim([0 1.5])
legend('\phi','h','\Sigma S_i + \Sigma S_{ij}','\Sigma S_{ij}')
xlabel('L [cm]'), ylabel('Sobol Indice'), title('Halton Sobol Indices vs. Length')

% % variance
% disp(VARy)
close(v1), close(v2), close(v3), 