%labb1
%% task3.1 q1
%ARMA: A(z)y_t=C(z)e_t
A1 = [ 1 -1.79 0.84];
C1 = [ 1 -0.18 -0.11];
arma1 = idpoly(A1, [], C1);
A2 = [ 1 -1.79];
C2 = [ 1 -0.18 -0.11];
arma2 = idpoly(A2, [], C2);
arma1.a;        %ger A(z) koeffs
%ger poles and zeros
pzmap(arma2)
subplot(121)
pzmap(arma1)
subplot(122)
pzmap(arma2)    %instabil

%creates noise with varaince sigma2
sigma2 = 1.5;
N = 200;
e = sqrt(sigma2) * randn(N, 1);
y_trend = filter(arma1.c, arma1.a, e);
y_per = filter(arma2.c, arma2.a, e);

figure
subplot(211)
plot(y_trend)
subplot(212)
plot(y_per)        %instabil
%% q2
m = 100; %cant calc higher than ~30
%computes covariance for arma, assumes sigma = 1
rtheo = kovarians(arma1.c, arma1.a, m);
figure
stem(0:m, rtheo * sigma2, 'b')
hold on
rest = covf(y_trend ,m+1);
stem(0:m, rest, 'r')
legend('Theoretical','with noice')
%% q3
my_plotter(y_trend)
na = 2;     nc = 2;
data = iddata(y_trend);
%LS based for ar
ar_model = arx(y_trend, [na]);
%error residual by inverted filter
ehatAR = filter(ar_model.a, ar_model.c, y_trend);
%ML based method for arma
arma_model = armax(y_trend, [na nc]);
%error residual by inverted filter
ehatARMA = filter(arma_model.a, arma_model.c, y_trend);
present(ar_model);
present(arma_model);
my_plotter(ehatAR)
my_plotter(ehatARMA)

%% 3.2 no q
n = 500;
A_star = [1 -1.35 0.43];
sigma2 = 4;
noise = sqrt(sigma2) * randn(n + 100, 1);
y = filter(1, A_star, noise);
y = y(101:end);     %corupted, not causal
subplot(211)
plot(y)
title('y')
subplot(212)
plot(noise)
title('noise')

nest = floor(2/3*n);             %floor avrundar ner till n?rmaste heltal
yest = iddata(y(1:nest));        %iddata ?r bara en klass?
yval = iddata(y(nest+1:end));
NN = [1: 10]';
V = arxstruc(yest, yval, NN);
n_order = selstruc(V, 0);
n_aic = selstruc(V, 'aic');

%% q4
n = 500;
sigma2 = 4;
Ar = [1 -1.35 0.43];
NN = [1: 10]';
for i= 1:100
    noise = sqrt(sigma2) * randn(n+100,1); 
    y = filter(1,Ar,noise);
    nest = floor(2/3*n);
    yest = iddata(y(1:nest)); 
    yval = iddata(y(nest+1:end));
    V = arxstruc(yest,yval,NN); 
    n_order(i) = selstruc(V,0);
    n_aic(i) = selstruc(V, 'aic');
end
subplot(211)
histogram(n_order)
title('LS')
subplot(212)
histogram(n_aic)
title('AIC')

armodel = arx(y, n_order(end));
armodel.NoiseVariance
armodel.CovarianceMatrix
present(armodel)

%% 3.3 q5
load data.dat

ar1_model = arx(data, 1);

figure
rar1 = resid(ar1_model, data);
present(ar1_model);

%% other orders
ar2_model = arx(data, 4);
ar3_model = arx(data, 5);
figure
rar2 = resid(ar2_model, data);
figure
rar3 = resid(ar3_model, data);
present(ar2_model);
present(ar3_model);

%% q6
am11_model = armax(data, [1 1]);
rar1 = resid(am11_model, data);
present(am11_model);

%% 3.4 q7
A_star = [1 -1.5 0.7];
C = [1 zeros(1, 11) -0.5];
A12 = [1 zeros(1, 11) -1];
A_star = conv(A_star, A12);
e = randn(600 , 1);
y = filter(C, A_star, e);
y = y(100 : end);
plot(y)
my_plotter(y, 100)

%% q8
y_s = filter(A12,1,y);
%data = iddata(y_s);
my_plotter(y_s)
%% q9
A = [1 -1 -1]
%A12 = [1 zeros(1, 11) -1];
%A_star = conv(A, A12);
B = [];
%C = [1 zeros(1, 11) -0.5];
%C = [1 -1 ]
%model_init = idpoly(A12,B,C);
model_init = idpoly( A, [], []);
model_armax = pem(y_s, model_init)
y_per = filter(model_armax.a, model_armax.c, y_s);
my_plotter(y_per)
%% q10
%data = iddata(y);
model_init = idpoly([1 0 0], [], [1 zeros(1, 12)]);
model_init.Structure.c.Free = [zeros(1, 12), 1];
model_armax = pem(y_s, model_init)
y_arma = filter(model_armax.a, model_armax.c, y_s);
my_plotter(y_arma)

%% 3.5
load('svedala.mat') 
plot(svedala)
my_plotter(svedala, 200)
%% removes trend
Adif = [1 -1];
y_trend = filter(Adif, 1, svedala);
y_trend = y_trend(10:end);
my_plotter(y_trend, 100)
%% removes periodisity
A24 = [1 zeros(1, 23) -1];
y_per = filter(A24, 1, y_trend);
y_per = y_per(50:end);
my_plotter(y_per, 100);
%present(y2);
%% models the process without trends and periods
A = [1 0];

model_init = idpoly(A, [], []);
%C = [1 zeros(1,23) 0];
%model_init.Structure.c.Free = [zeros(1, 24) 1];
model_armax = pem(y_per, model_init);
present(model_armax)
y_arma = filter(model_armax.a, model_armax.c, y_per);
y_arma = y_arma(100:end);
my_plotter(y_arma, 100);
%%
figure
whitenessTest(y_arma)
figure 
plot(y_arma)
