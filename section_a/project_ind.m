function [T] = project_ind(xx,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project
% Author: Danial Faghihi, University at Buffalo
% Edited: Alden Noel, UB Grad         
% INPUTS:
% xx = [x1, x2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = xx(1);
x2 = xx(2);
x3 = x;
% x3 = xx(3);
% x3 = linspace(10,66,N)';
% x3 = x3(N);
a = 0.95;     % cm
b = 0.95;     % cm
L = 70.0;     % cm
k = 2.37;     % W/cm C
T_amb = 21.29; % C
gamma = sqrt((2*(a+b)*x2)/(a*b*k));
c1 = -(x1/(k*gamma))*((exp(gamma*L)*(x2+k*gamma))/(exp(-gamma*L)*(x2-k*gamma)+exp(gamma*L)*(x2+k*gamma)));
c2 = x1/(k*gamma) + c1;

term1 = c1*exp(-gamma*x3);
term2 = c2*exp(gamma*x3);
term3 = T_amb;

T = term1 + term2 + term3;


end
