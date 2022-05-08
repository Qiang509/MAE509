close all;
clear all;
clc;

a=.95;
b=.95;
L=70;
h= .00191;
phi=-18.41;
k=2.37;
Tamb=21.29;

gamma=sqrt((2*(a+b)*h)/(a*b*k));
c1=(-phi/(k*gamma))*(((exp(gamma*L))*(h+(k*gamma)))/  (exp(-gamma*L)*(h-(k*gamma))+(exp(gamma*L)*(h+(k*gamma)))));
c2=(phi/(k*gamma))+c1;

% x=[0 10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
%x=linspace(10, 66, 1000);
udata = [96.1 80.12 67.66 57.96 50.90 44.84 39.75 36.16 33.31 31.15 29.28 27.88 27.18 26.40 25.86];
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  STD_variables=[0.2 0.5 0.8 .45 .32 .15 .7 .65 .54 .48 .84 .56 .74 .36 .75];
  error=STD_variables.*2.5;
x=linspace(10, 66,length(udata));

T=(c1*exp(-gamma.*x))+((c2*exp(gamma.*x)))+Tamb;
figure()
plot(x,T,'-k.','MarkerSize',10);
hold on
errorbar(xdata,udata,error,'--r')
grid on
grid minor 
title('Temperature distribution on initial bar')
ylabel('Temoerature [degrees C]')
xlabel('Distance (cm)')
legend('Model Results', 'Given Data')
hold off


