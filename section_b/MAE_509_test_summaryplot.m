%MAE_509_test_summaryplot.m


load T1.mat q1guess q2guess chain results
q1guess
q2guess
chain1 = chain;
T_1 = results;

load T2.mat q1guess q2guess chain results
q1guess
q2guess
chain2 = chain;
T_2 = results;

load T3.mat q1guess q2guess chain results
chain3 = chain;
q1guess
q2guess
T_3 = results;

load T4.mat q1guess q2guess chain results
chain4 = chain;
q1guess
q2guess
T_4 = results;

load T5.mat q1guess q2guess chain results
chain5 = chain;
q1guess
q2guess
T_5 = results;

figure(1)
mcmcplot(chain1,[],T_1,'chainpanel');
%title('t1','t2')
%legend('n')
figure(2)
mcmcplot(chain2,[],T_2,'chainpanel');
%title('t2')
figure(3)
mcmcplot(chain3,[],T_3,'chainpanel');
%title('t3')
figure(4)
mcmcplot(chain4,[],T_4,'chainpanel');
%title('t4')
figure(5)
mcmcplot(chain5,[],T_5,'chainpanel');
%title('t5')

a = mean(chain1((1000:10000),1));
b = mean(chain2((1000:10000),1));
c = mean(chain3((1000:10000),1));
d = mean(chain4((1000:10000),1));
e = mean(chain5((1000:10000),1));

Q1 = (a+b+c+d+e)/5

aa = mean(chain1((1000:10000),2));
bb = mean(chain2((1000:10000),2));
cc = mean(chain3((1000:10000),2));
dd = mean(chain4((1000:10000),2));
ee = mean(chain5((1000:10000),2));

h1 = (aa+bb+cc+dd+ee)/5
