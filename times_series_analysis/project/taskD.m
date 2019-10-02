%modelling for Kassala
% Apply the better of your models to the data from Kassala. 
% Are the prediction errors similar? Are the results improved 
% if you re-estimate the model parameters based on the different 
% measurements? Locate the two measurement stations on a map and 
% discuss some possible reasons for the data to behave differently
% at the two locations.
load proj18
nvdi = Kassala.nvdi;
nvdi = double(nvdi/128) - 1;
nvdi = log(nvdi);%log was most normal distributed in jbtest and lillietest
MEAN_KASSALA  = mean(nvdi);
nvdi = nvdi - MEAN_KASSALA;
model_cut = floor(length(nvdi)*0.7);
m_Kassala = nvdi(1:model_cut);
v_Kassala = nvdi(model_cut + 1:end);

clear nvdi
%% PLOT INITIAL DATA
plot(m_Kassala)
my_plotter(m_Kassala)

%% Use A and C ploynomials form ELGeneina
A36 = [1 zeros(1,35) -1];
A_old = conv([1, -0.6834], A36);
C_old = 1;
Kassala_veg = filter(A_old, C_old, m_Kassala);
Kassala_veg = Kassala_veg(100:end);
my_plotter(Kassala_veg);
% we get white noise
figure
plot(Kassala_veg)
figure
whitenessTest(Kassala_veg)
clear Kassala_veg
%% MODELING THE VEGITATION
% removes sesonality
A36 = [1 zeros(1, 35) -1];
m_Kass = filter(A36, 1, m_Kassala);
m_Kass = m_Kass(50:end);
my_plotter(m_Kass);
%% get new A and C
A = [1 0]; %1 eller 2?
%C = [1 zeros(1,35) 0];
C = 1;
kas_mod_init = idpoly(A, [], C);
%kas_mod_init.Structure.c.Free = [zeros(1, 36) 1];
kas_model = pem(m_Kass, kas_mod_init);
present(kas_model)
SARIMA_kas = filter(kas_model.a, kas_model.c, m_Kass);
SARIMA_kas = SARIMA_kas(100:end);
my_plotter(SARIMA_kas);
% we get white noise
figure
plot(SARIMA_kas)
figure
whitenessTest(SARIMA_kas)
clear A C kas_mod_init SARIMA_kas
%% Preds for old model HOW do i judge the kvalle of a pred?
y = v_Kassala;
k = 5;
%pred using the ELGeneina model
A = A_old;%conv(vegitation_model.a, A36);
C = C_old;%vegitation_model.c;
[CS,AS] = equalLength(C, A);
[Fk,Gk1] = deconv(conv([1,zeros(1,k-1)],CS),AS);
pred_old = filter(Gk1,C,y);
%pred using new model
A = conv(kas_model.a, A36);
C = kas_model.c;
[CS,AS] = equalLength(C, A);
[Fk,Gk2] = deconv(conv([1,zeros(1,k-1)],CS),AS);
pred_new = filter(Gk2,C,y);

pred_old = pred_old(length(Gk1)+5: end);
pred_new = pred_new(length(Gk1)+5: end);

naivepred = zeros(1,length(y));
for i = 1:length(y)-k
    naivepred(i+k) = y(i);
end
y = y(length(Gk1)+5: end);
naivepred = naivepred(length(Gk1)+5: end);

%take away the transforms
pred_old = exp(pred_old+MEAN_KASSALA);
pred_new = exp(pred_new+MEAN_KASSALA);
y = exp(y+MEAN_KASSALA);
naivepred = exp(naivepred+MEAN_KASSALA);


%% plots
%plotta preds + ACF/PACF f?r felen samt ber?kna felvariansen
figure
subplot(2,1,1)
title('The 5-step predictions')

hold on
plot(pred_old,'r--')
plot(y,'b')
plot(naivepred,'g.-')
axis tight
legend('El-Genieinas model', 'validation data', 'naive prediction')
hold off
subplot(2,1,2)
hold on
plot(pred_new,'r--')
plot(y,'b')
plot(naivepred,'g.-')
axis tight
legend('Kassalas model', 'validation data', 'naive prediction')
hold off

errornaive = y-naivepred';
errorprediction1 = y - pred_old;
errorprediction2 = y - pred_new;
figure
subplot(3,1,1)
acf(errorprediction1,50,0.05,1);
title('ACF for the error in the 5-step prediction with El-Genieinas model')
subplot(3,1,2)
acf(errorprediction2,50,0.05,1);
title('ACF for the error in the 5-step prediction with Kassalas model')
subplot(3,1,3)
acf(errornaive,50,0.05,1);
title('ACF for the error in the 5-step prediction with the naive prediction')
%0.0705
variances  = [var(errorprediction1), var(errorprediction2), var(errornaive)]
