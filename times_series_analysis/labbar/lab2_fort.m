%3.3
load svedala
y = svedala;
A= [1 -1.79 0.84]; 
C= [1 -0.18 -0.11];
k = 1;
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhatk = filter(Gk,C,y);
yhatk = yhatk(101:end);

totaldata = [y(101:end) yhatk];
figure
plot(totaldata)


for i = 101:length(y)
    error(i) = y(i)-yhatk(i-100);
end

figure
plot(error)
%figure
%whitenessTest(error)

my_plotter(error)
mean(error)

%% 3.4
load sturup
tempsturup = sturup;
%figure
%plot(tempsturup)
Au = [ 1 -1.49 0.57];
Bu = [ 0 0 0 0.28 -0.26 ]; %delay d is 3 samples
C = [1];
Gyk = Gk;
Fyk = Fk;
k = 3;
BFk = conv(Bu,Fyk);
%svedala pred
y_hat = filter(BFk, C, sturup) + filter(Gyk, C, svedala);
y_hat = y_hat(101:end);

[BFs,Cu] = equalLength(BFk,C);
[FKu,GKu] = deconv(conv([1,zeros(1,k-1)],BFs),Cu);

u_h = filter(GKu, C, tempsturup);
y1_h = filter(Gyk, C, svedala);

y_h = u_h + y1_h;
y_h=y_h(101:end);

plot(svedala(101:end));
hold on
plot(y_hat)

error2 = y_h - svedala(101:end);
variance2 = std(error2).^2;
figure
whitenessTest(error2)
my_plotter(error2)

m_err2 = mean(error2)
Expectation2 = sqrt(variance2) + m_err2;
figure
plot(error2)


%% 3.5
load svedala
plot(svedala)
my_plotter(svedala, 100)
% removes trend
% Adif = [1 -1];
% y_trend = filter(Adif, 1, svedala);
% y_trend = y_trend(10:end);
% my_plotter(y_trend, 100)
% removes periodisity
A24 = [1 zeros(1, 23) -1];
y_per = filter(A24, 1, svedala);
y_per = y_per(50:end);
my_plotter(y_per, 100);
%present(y2);
%% models the process without trends and periods
A = [1 0 0];
S = 24;
C = [1 zeros(1,23) 0];
model_init = idpoly(A, [], C);
model_init.Structure.c.Free = [1 zeros(1, 23) 1];
model_sarima = pem(y_per, model_init);
present(model_sarima)
y_sarima = filter(model_sarima.a, model_sarima.c, y_per);
y_sarima = y_sarima(100:end);
my_plotter(y_sarima, 100);

figure
whitenessTest(y_sarima)
figure 
plot(y_sarima)
%% kombines the filters to predict svedala
y = svedala;

A24 = [1 zeros(1, 23) -1];
A1 = model_sarima.a;
A = conv(A1, A24);
C = model_sarima.c;
k = 1;
[CS,AS] = equalLength(C, A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
yhatk = filter(Gk,C,y);

plot(y);
hold on
plot(yhatk)

err = yhatk - y;
err = err(100:end);
variance = std(err).^2

figure
whitenessTest(err)
my_plotter(err)
%%
m_err = mean(abs(err))
Expectation2 = sqrt(variance) + m_err
figure
plot(err)

