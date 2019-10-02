%% 3.4 q6
load svedala94
y = svedala94;
T = linspace(datenum(1994, 1, 1), datenum(1994,12,31), length(y));
plot(T,y);
datetick('x');
my_plotter(y)
%% remove trend
Adif = [1 0 0 0 0 0 -1];
y_diff = filter(Adif, 1, y);
y_diff = y_diff(10:end);
plot(T(10:end), y_diff)
datetick('x');
%my_plotter(y_diff, 100)
%%
% Use armax to estimate an ARMA(2,2)-process for the differentiated data.
% Do this (i) for the entire year, (ii) for January-March, and (iii) for
% June-August. Compare the different estimated parameters.

%index tokigt n?r vi kastat bort datapunkter?
th = armax(y_diff,[2 2]);
th_winter = armax(y_diff(1:540),[2 2]);
th_summer = armax(y_diff(907:1458),[2 2]);
%% summer is not to white
arma = filter(th.a, th.c, y_diff);
arma_winter = filter(th_winter.a, th_winter.c, y_diff(1:540));
arma_summer = filter(th_summer.a, th_summer.c, y_diff(907:1458));

%%
plot(arma_summer)
my_plotter(arma_summer)
%% looking at the estimated parameters
th.a, th.c
th_summer.a, th_summer.c
th_winter.a, th_winter.c

%% q7
th0 = [th_winter.A(2:end) th_winter.C(2:end)];
%th0 = [0 0 0 0];
lambda = 0.99;  %0.99 best
[thr, yhat] = rarmax(y_diff,[2 2],'ff',lambda,th0);
figure
subplot(311)
plot(T,svedala94);
datetick('x')
subplot(312)
plot(thr(:,1:2))
hold on
plot(repmat(th_winter.A(2:end),[length(thr) 1]),'g:');
plot(repmat(th_summer.A(2:end),[length(thr) 1]),'y:');
axis tight
hold off

subplot(313)
plot(thr(:,3:end))
hold on
plot(repmat(th_winter.C(2:end),[length(thr) 1]),'g:');
plot(repmat(th_summer.C(2:end),[length(thr) 1]),'y:');
axis tight
hold off
%% 3.5 q8
load svedala94
y = svedala94(850:1100);

t = (1:length(y))';
U = [sin(2*pi*t/6) cos(2*pi*t/6)];
Z = iddata(y,U);
model = [3 [1 1] 4 [0 0]];
%[na[nb1nb2]nc[nk1nk2]];
thx = armax(Z,model);
c2mt = cell2mat(thx.b)';
%f?r b?ttre plot
y = y - mean(y);
plot(y)
hold on
plot(U * c2mt, 'g')

%% q9 
U=[sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z=iddata(y,U);
m0=cell2mat([thx.A(2:end) thx.B 0 thx.C(2:end)]);
Re = 1*diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re, m0);
% vad ?r detta
sved  = svedala94(861:1100);
plot(yhat)
% figure 
% plot(sved)
%%
%motsvarar var ettan ?r i Re, rimligt? ger helt annat svar
m = thr(:, 6); 
a = thr(end, 4);
b = thr(end, 5);
y_mean = m + a*U(:, 1) + b*U(:,2);
y_mean = [0; y_mean(1: end - 1)];
%vad ?r detta?
%figure
hold on
plot(y_mean)


%% q10 g?r om fast med y = y - y(1);
load svedala94
%varf?r minskas amplituden p? U*thx.b n?r intervaller f?rl?ngs?
y = svedala94;
y = y - y(1);
t = (1:length(y))';
U = [sin(2*pi*t/6) cos(2*pi*t/6)];
Z = iddata(y,U);
model = [3 [1 1] 4 [0 0]];
%[na[nb1nb2]nc[nk1nk2]];
thx = armax(Z,model);
c2mt = cell2mat(thx.b)';

y_diff = filter([1 -1], 1, y);
y_diff = y_diff(11:end);
plot(y_diff)
hold on
plot(U(11:end,:) * c2mt, 'g')
axis tight
%% 
U=[sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z=iddata(y,U);
m0=cell2mat([thx.A(2:end) thx.B 0 thx.C(2:end)]);
%%
Re = diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re, m0);
plot(yhat)%similar to svedala94
err = svedala94-yhat;   %not white

%%
%motsvarar var ettan ?r i Re, rimligt? ger helt annat svar
m = thr(:, 6);  
a = thr(end, 4);
b = thr(end, 5);
y_mean = m + a*U(:, 1) + b*U(:,2);
y_mean = [0; y_mean(1: end - 1)];
%vad ?r detta?
plot(y_mean)

figure
plot(svedala94)


