% As an optional task, consider making the better of your models 
% recursive, such that the parameters of the model are allowed 
% to vary over time. Does this improve the prediction ability as 
% compared to your non-recursive model? When forming predictions 
% using the precipitation data as an external input, you may treat 
% future values of this data as known. Thus, even if needed,
% you do not need to predict the precipitation data 
% (but you are of course welcome to do so if you wish).
load proj18

nvdi = ElGeneina.nvdi;
nvdi = double(nvdi/128) - 1;
nvdi = log(nvdi);%log was most normal distributed in jbtest and lillietest
nvdi = nvdi - mean(nvdi);
% The good model I choose is the Regular one which had the order
% Arma(2,1)
A36 = [1 zeros(1, 35) -1];
nvdiP = filter(A36, 1, nvdi);
nvdiP = nvdiP(50:end);
my_plotter(nvdiP);

%% AR(1) after season is removed
na = 1;
model = [na];%,nb,nk]
n = 100;
lambdaline = linspace(.97, 1, n);
ls2 = zeros(n, 1);
for i = 1:n
    [Aest, yhat, CovAest, trash] = rarx(nvdiP, model, 'ff', lambdaline(i));
    ls2(i) = sum((nvdiP - yhat).^2);
end
plot(lambdaline, ls2) 
hold on
axis tight
argmin = find((ls2-min(ls2))<=0)
opt_lambda = lambdaline(argmin)
plot(opt_lambda, ls2(argmin),'*')
%%
guess = -0.6784;
lambda = opt_lambda; %0.95 ives a good estimate with reasonable variance
%use recursiveARX instead of rarx on newer matlab
[Aest, yhat, covAest, yprev] = rarx(nvdiP, model, 'ff', lambda);

figure
subplot(2,1,1)
plot(yhat);
subplot(2,1,2)
plot(Aest, 'b');
hold on
plot([0 600], [guess guess], 'r-.');
%% PREDICTS THE VEGITATION
% The ARMA-prediction is not great
y = nvdi(end-199:end);
A = [ones(200, 1), Aest(end-199:end)];
C = 1;
k = 2;
A1 = conv(A(1,:), A36);
[CS,AS] = equalLength(C, A1);
FK = [];
GK = [];
RAR_pred = [];
%Funkar inte, m?ste h?rleda ekvationerna f?r uppdatering av G
for i = 1:10:length(y)-k
    Ai = conv(A(i,:), A36);
    [Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),Ai);
    RAR_pred = [RAR_pred; filter(Gk,C,y(i:i+k-1))];
end

figure
plot(y, '-.')
hold on
plot(RAR_pred, '*')
legend('nvdi', 'recursive prediction')

%%
ARMA_pred = filter(Gk,C,y);

plot(y(50:end), 'r-.')
hold on
plot(ARMA_pred(50:end))

%% PREDICTS THE VEGITATION
% The ARMA-prediction is not great
y = nvdiP(end-200:end);
A = Aest(end-200:end);
C = 1;
k = 5;
RAR_pred = [];
for i = 1:k:length(y)-k
    Ai = A(i);
    [CS,AS] = equalLength(C, Ai);
    [Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
    RAR_pred = [RAR_pred; filter(Gk,C,y(i:i+k-1))];
end
figure
plot(y, 'r-.')
hold on
plot(RAR_pred)

err = y(50:end) - ARMA_pred(50:end);
variance = var(err);
figure
whitenessTest(err)
my_plotter(err)
figure
plot(err)
m_err = mean(err)
Expectation = sqrt(variance) + m_err
