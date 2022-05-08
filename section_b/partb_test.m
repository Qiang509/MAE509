  clc
  clear all;
  close all;

  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  u_amb = 21.29;
  stD = [0.2,0.5,0.8,0.45,0.32,0.15,0.7,0.65,0.54,0.48,0.84,0.56,0.74,0.36,0.75]*2.5;
  a = 0.95;   % cm
  b = 0.95;   % cm
  L = 70.0;   % cm
  k = 2.37;   % W/cm C
  h = 0.00191;
  Q = -18.41; 
  %Q = -20.41
  %h =0.0191;
  n = 15;
  p = 2;
  
  gamma = sqrt(2*(a+b)*h/(a*b*k));
  gamma_h = (1/(2*h))*gamma;
  f1 = exp(gamma*L)*(h + k*gamma);
  f2 = exp(-gamma*L)*(h - k*gamma);
  f3 = f1/(f2 + f1);
  f1_h = exp(gamma*L)*(gamma_h*L*(h+k*gamma) + 1 + k*gamma_h);
  f2_h = exp(-gamma*L)*(-gamma_h*L*(h-k*gamma) + 1 - k*gamma_h);
  c1 = -Q*f3/(k*gamma);
  c2 = Q/(k*gamma) + c1;
  f4 = Q/(k*gamma*gamma);
  den2 = (f1+f2)^2;
  f3_h = (f1_h*(f1+f2) - f1*(f1_h+f2_h))/den2;
  c1_h = f4*gamma_h*f3 - (Q/(k*gamma))*f3_h;
  c2_h = -f4*gamma_h + c1_h;
  c1_Q = -(1/(k*gamma))*f3;
  c2_Q = (1/(k*gamma)) + c1_Q;

  uvals_data = c1*exp(-gamma*xdata) + c2*exp(gamma*xdata) + u_amb;
  uvals_Q_data = c1_Q*exp(-gamma*xdata) + c2_Q*exp(gamma*xdata);
  uvals_h_data = c1_h*exp(-gamma*xdata) + c2_h*exp(gamma*xdata) + gamma_h*xdata.*(-c1*exp(-gamma*xdata) + c2*exp(gamma*xdata));
  
  plot(xdata,uvals_data)
  hold on
  errorbar(xdata,uvals_data,stD,'or','LineWidth',2)
  hold off
