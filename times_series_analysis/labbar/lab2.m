% tidserie lab2
% 3.1
n = 500;
A1 = [1 -.65]; A2 = [1 .90 .78];
C = 1 ; B = [0 0 0 0 .4];
e = sqrt(1.5) * randn(n + 100, 1);
w = sqrt(2) * randn(n + 200, 1);
A3 = [1 .5]; C3 = [1 -.3 .2];
u = filter(C3 , A3, w); u = u(101: end);
y = filter(C, A1, e) + filter(B, A2, u);
u = u(101:end); y = y(101: end);

clear A1 A2 C B e  A3 C3
%plot(u)
%figure
%plot(y)
%% q1 ut = C/A * upw pre whitening
 plot(u)
my_plotter(u)
%% a simple AR(1) passed the whiteness test
A = [1 0];
model_init = idpoly(A, [], []);
model_ar = pem(u, model_init);
present(model_ar)
w = filter(model_ar.a, model_ar.c, u);
w = w(100:end);
my_plotter(w, 100);
figure
whitenessTest(w)
figure 
plot(w)

%% q2 pw p? y
ypw = filter(model_ar.a, model_ar.c, y);
ypw = ypw(100:end);

M = 40; 
stem(-M:M, crosscorr(w , ypw ,M));
title('Cross-correlation-function'), xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
%% q2 fort calc H. v is not white but almost uncorrelated
A2 = [1 0 0];
B = [0 0 0 0 1];
Mi = idpoly([1], [B], [], [], [A2]);
Mi.Structure.b.Free = [0 0 0 0 1];
zpw = iddata(ypw, w);

Mba2 = pem(zpw, Mi); 
present(Mba2)
vhat = resid(Mba2, zpw,'corr')
v = vhat.OutputData;
stem(-M:M, crosscorr(w , v ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
my_plotter(v)
figure
whitenessTest(v)
%% q3  tittar p? x
x = y - filter(Mba2.b, Mba2.f, u);
x = x(100:end);
my_plotter(x)
%% x can be modeld fine by AR1 and is uncorr with u
A1 = [1 0];
model_init = idpoly(A1, [], []);
model_ar = pem(x, model_init);
present(model_ar)
x_ar = filter(model_ar.a, 1, x);
x_ar = x_ar(100:end);
my_plotter(x_ar, 100);
figure
whitenessTest(x_ar)

figure
stem(-M:M, crosscorr(u , x ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
%% q4
A1 = model_ar.a;
A2 = Mba2.f;
B = Mba2.b;
C = 1;  %model_ar.c
M_i = idpoly(1 ,B, C, A1, A2);
z = iddata(y, u);
MboxJ = pem(z, M_i);
present(MboxJ)
ehat = resid(MboxJ, z);

e = ehat.OutputData;
figure
stem(-M:M, crosscorr(u , e ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
my_plotter(e)
figure
whitenessTest(e)

%% 3.2 q5 hairdryer

load('tork.dat')
tork = tork - repmat(mean(tork), length(tork), 1);
u = tork(:, 2);%input
y = tork(:, 1);%output 
z = iddata(y,u);
plot(z(1:300))
%% ut = C/A * w  vi kan pw u med ar1
my_plotter(u)
A = [1 0];
model_init = idpoly(A, [], []);
model_ar = pem(u, model_init);
present(model_ar)
w = filter(model_ar.a, model_ar.c, u);
w = w(100:end);
my_plotter(w, 100);
figure
whitenessTest(w)
% figure 
% plot(w)
%% kollar p? korrelationen mellan w och epsilon
epsilon = filter(model_ar.a, model_ar.c, y);
epsilon = epsilon(100:end);

M = 40; 
stem(-M:M, crosscorr(w, epsilon, M)); %p? vilker h?ll ska den st?
title('Cross-correlation-function'), xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
%% tar fram H med f?rra cr_cor. inte white men n?stan uncorr
A2 = [1 0 0];
B = [0 0 0 1 1 1];
Mi = idpoly([1], [B], [], [], [A2]);
Mi.Structure.b.Free = [0 0 0 1 1 1 ];
zpw = iddata(epsilon, w);

Mba2 = pem(zpw, Mi); 
present(Mba2)
vhat = resid(Mba2, zpw,'corr')
v = vhat.OutputData;
figure
stem(-M:M, crosscorr(w , v ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
my_plotter(v)
% figure
% whitenessTest(v)
%%  Beh?ver vi forts?tta?
x = y - filter(Mba2.b, Mba2.f, u);
x = x(100:end);
my_plotter(x)
%% vitt men korrelerade
A1 = [1 0];
model_init = idpoly(A1, [], []);
model_ar = pem(x, model_init);
present(model_ar)
x_ar = filter(model_ar.a, 1, x);
x_ar = x_ar(100:end);
my_plotter(x_ar, 100);
figure
whitenessTest(x_ar)

figure
stem(-M:M, crosscorr(u , x ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
%% kanske b?ttre om vi b?rjar med ex.[1 0 0] polynom
A1 = model_ar.a;
A2 = Mba2.f;
B = Mba2.b;
C = 1;  %model_ar.c
M_i = idpoly(1 ,B, C, A1, A2);
z = iddata(y, u);
MboxJ = pem(z, M_i);
present(MboxJ)
ehat = resid(MboxJ, z);

e = ehat.OutputData;
figure
stem(-M:M, crosscorr(u , e ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
my_plotter(e)
figure
whitenessTest(e)






