load proj18
% El-Geneina, standard data modifiering
%vegitation . Ta bort outliers??
nvdi = ElGeneina.nvdi;
nvdi = double(nvdi/127.5) - 1;
nvdi = log(nvdi);%log was most normal distributed in jbtest and lillietest
MEAN_VALUE_NVDI = mean(nvdi);
nvdi = nvdi - MEAN_VALUE_NVDI;
model_cut = floor(length(nvdi)*0.7);
m_data = nvdi(1:model_cut);
v_data = nvdi(model_cut + 1:end);
%% the rain
cut_off = length(zsave) - length(nvdi);
rain = zsave(cut_off + 1:end)';
rain = rain - mean(rain);
m_rain = rain(1:model_cut);
v_rain = rain(model_cut + 1:end);

clear nvdi rain model_cut cut_off
%% PLOT INITIAL DATA
plot(iddata(m_data, m_rain))
my_plotter(m_data,100)
my_plotter(m_rain,100)
%% MODELING THE VEGITATION
% removes sesonality
A36 = [1 zeros(1, 35) -1];
m_dataP = filter(A36, 1, m_data);
m_dataP = m_dataP(50:end);
my_plotter(m_dataP);
%% Estimates A and C ploynomials
A = [1 0]; %1 eller 2?
C = [1 zeros(1,35) 0]; %b?ttre utan, lol
%C = 1;
veg_mod_init = idpoly(A, [], C);%
veg_mod_init.Structure.c.Free = [zeros(1, 36) 1];
vegitation_model = pem(m_dataP, veg_mod_init);
present(vegitation_model)
SARIMA_veg = filter(vegitation_model.a, vegitation_model.c, m_dataP);
SARIMA_veg = SARIMA_veg(100:end);
my_plotter(SARIMA_veg);
% we get white noise
figure
plot(SARIMA_veg)
figure
whitenessTest(SARIMA_veg)
clear A C veg_mod_init SARIMA_veg 
%% PREDICTS THE VEGITATION
y = v_data;
k = 5;
% Model 1
A = conv([1, -0.6834], A36);
C = 1;
[CS,AS] = equalLength(C, A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
ARMA_pred1 = filter(Gk,C,y);
ARMA_pred1 = ARMA_pred1(50: end);
% model 2
A = conv([1, -.6618], A36);
C = [1, zeros(1,35), -0.7348];
[CS,AS] = equalLength(C, A);
[Fk,Gk] = deconv(conv([1,zeros(1,k-1)],CS),AS);
ARMA_pred2 = filter(Gk,C,y);
ARMA_pred2 = ARMA_pred2(50: end);

naivepred = zeros(1,length(y));
for i = 1:length(y)-k
    naivepred(i+k) = y(i);
end
y = y(50: end);
naivepred = naivepred(50: end);

%take away the transforms
ARMA_pred1 = exp(ARMA_pred1+MEAN_VALUE_NVDI);
ARMA_pred2 = exp(ARMA_pred2+MEAN_VALUE_NVDI);
y = exp(y+MEAN_VALUE_NVDI);
naivepred = exp(naivepred+MEAN_VALUE_NVDI);


%% plots
%plotta preds + ACF/PACF f?r felen samt ber?kna felvariansen
figure
subplot(2,1,1)
title('The 5-step predictions')

hold on
plot(ARMA_pred1,'r--')
plot(y,'b')
plot(naivepred,'g.-')
axis tight
legend('Model1s predicion', 'validation data', 'naive prediction')
hold off
subplot(2,1,2)
hold on
plot(ARMA_pred2,'r--')
plot(y,'b')
plot(naivepred,'g.-')
axis tight
legend('Model2s predicion', 'validation data', 'naive prediction')
hold off

errornaive = y-naivepred';
errorprediction1 = y - ARMA_pred1;
errorprediction2 = y - ARMA_pred2;
figure
subplot(3,1,1)
acf(errorprediction1,50,0.05,1);
title('ACF for the error in the 5-step prediction with model 1')
subplot(3,1,2)
acf(errorprediction2,50,0.05,1);
title('ACF for the error in the 5-step prediction with model 2')
subplot(3,1,3)
acf(errornaive,50,0.05,1);
title('ACF for the error in the 5-step prediction with the naive prediction')
%0.0705
variances  = [var(errorprediction1), var(errorprediction2), var(errornaive)]
% figure
% whitenessTest(err)
% my_plotter(errorprediction)
% my_plotter(errornaive)
% figure
% plot(err)
% m_err = mean(err)
% Expectation = sqrt(variance_pred) + m_err
%clear y  C CS AS err m_err Expectation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARMAX: MODEL THE VEGITATION WITH RAIN AS INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREWHITE THE RAIN
%Removing the season
A36 = [1 zeros(1, 35) -0.6];
m_rp = filter(A36, 1, m_rain);
m_rp = m_rp(50:end);
my_plotter(m_rp);

%% Modeling the rain
%this model works as a treat, somehow
A = [1 0 0];
C = [1 zeros(1,36)];
pw_mod_i = idpoly(A, [], C);
pw_mod_i.Structure.c.Free = [zeros(1, 36) 1];
pw_model = pem(m_rp, pw_mod_i);
present(pw_model)
w = filter(pw_model.a, pw_model.c, m_rp);
w = w(100:end);
A3 = conv(pw_model.a, A36);
C3 = pw_model.c;
my_plotter(w, 100);
figure
whitenessTest(w)
figure 
plot(w)
clear A C pw_mod_i pw_model
%% FIND THE TRANSFER FUNCTION, tittar vi innanf?r CI syns d=3, r=2, s=2
epsilon = filter(A3, C3, m_data);
cut = length(epsilon) - length(w);
epsilon = epsilon(cut + 1:end);
M = 40; 
stem(-M:M, crosscorr(w, epsilon, M)); 
title('Cross-correlation-function'), xlabel('Lag')
hold on
%r?tt conf int%
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
clear cut
%% H blir d? bara 1?
A2 = [1];%[1 0 0];
B = [0 0 1];%[0 0 0 1 1 1];
H_init = idpoly([1], [B], [], [], [A2]);
%Mi.Structure.b.Free = [0 0 0 1 1 1 ];
zpw = iddata(m_data, m_rain);
H_model = pem(zpw, H_init); 
present(H_model)
vhat = resid(H_model, zpw,'corr');
%viktiga polynom
A2  = H_model.f;
Bzd = H_model.b;
v_hat = vhat.OutputData;
figure
stem(-M:M, crosscorr(w, v_hat ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
my_plotter(v_hat)
% figure
% whitenessTest(v)
clear H_init zpw H_model vhat v_hat
%% MODELING x = C1/A1 et = vegitation - B/A2 rain 
x = m_data - filter(Bzd, A2, m_rain);
x = x(100:end);
my_plotter(x)
%% remove sesonality
A36 = [1 zeros(1, 35) -.5];
x_s = filter(A36, 1, x);
x_s = x_s(50:end);
my_plotter(x_s);

%% MAKE THE ARMA MODEL AND OBTAIN A1, C1
A1 = [1 0];
C1 = [1];
model_init = idpoly(A1, [], C1);
AR_model = pem(x_s, model_init);
present(AR_model)
x_ar = filter(AR_model.a, AR_model.c, x);
x_ar = x_ar(100:end);
A1 = conv(AR_model.a, A36);
C1 = AR_model.c;
my_plotter(x_ar, 100);
figure
whitenessTest(x_ar)
figure
stem(-M:M, crosscorr(m_rain ,x ,M));
hold on
plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1) , '--' )
plot(-M:M, -2/sqrt(length(w))*ones(1, 2*M+1), '--')
hold off
clear model_init AR_model x_ar
%% PUT THE MODEL TOGETHER
%ange polynomen bin?rt h?r, beh?ll bara storleken
B = Bzd;
BJ_init = idpoly(1 ,B, C1, A1, A2);
BJ_init.Structure.d.Free = [1 1 zeros(1, 34) 1 1];
BJ_init.Structure.b.Free = [0 0 1];
z = iddata(m_data, m_rain);
MboxJ = pem(z, BJ_init);
present(MboxJ)
ehat = resid(MboxJ, z);
e_hat = ehat.OutputData;
figure
stem(-M:M, crosscorr(m_rain , e_hat ,M));
hold on

plot(-M:M, 2/sqrt(length(w)) * ones(1, 2*M+1), '--' )
plot(-M:M, -2/sqrt(length(w)) * ones(1, 2*M+1), '--')
hold off
my_plotter(e_hat)
figure
whitenessTest(e_hat)
clear BJ_init z MboxJ ehat e_hat
%% kombines the filters to predict
y = v_data;
u = v_rain;
%figure
%plot(tempsturup)
Au = conv(A1, A2);
Bu = conv(A1, Bzd); %fick d=s=r=0
C  = conv(C1, A2);
%G och F fr?n sarima. VIKTIGT ATT K ?R R?TT!!
Gyk = Gk;
Fyk = Fk;
% k = 1; we keep the choise from the arma prediction
BFk = conv(Bu, Fyk);
%ARMAX pred
% naive_pred = filter(BFk, C, u) + filter(Gyk, C, y);
% naive_pred = naive_pred(50:end);

[BFs, Cu]  = equalLength(BFk, C);
[FKu, GKu] = deconv(conv([1, zeros(1,k-1)], BFs), Cu);

u_h  = filter(GKu, C, u);
y1_h = filter(Gyk, C, y);

ARMAX_pred = y1_h; % + u_h;
ARMAX_pred = ARMAX_pred(50:end);

plot(y(50:end))
hold on
plot(ARMAX_pred)
title('ARMAX prediction')
%%
% figure
% plot(y(50:end));
% hold on
% plot(naive_pred)
% title('na?ve prediction')
figure
plot(y(50:end));
hold on
plot(ARMAX_pred)
title('ARMAX prediction')

err = y(50:end) - ARMAX_pred;
variance_pred = var(err);
figure
whitenessTest(err)
my_plotter(err)
figure
plot(err)
m_err = mean(err)
Expectation = sqrt(variance_pred) + m_err

clear u y Au Bu C Gyk Fyk BFk BFs Cu FKu GKu
clear u_h y1_h