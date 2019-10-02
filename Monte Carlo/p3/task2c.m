%% c
load('atlanticdata.mat')
F_inv = @(u,mu,beta) mu-beta.*(log(-log(u)));
T = 4200;
l = length(atlantic);
n = 20000;
[beta mu] = est_gumbel(atlantic);
tau_hat = F_inv(1-1/T,mu,beta);

t_th = zeros(n,1);
t_new = zeros(n,1);

for i = 1:n
    P0_hat = F_inv(rand(l,1),mu,beta);
    [x y] = est_gumbel(P0_hat);
    t_th(i) = F_inv(1-1/T,y,x);
    t_new(i) = F_inv(1-1/T,y,x)-tau_hat;
end

delta = mean(t_new);
t_th = t_th-delta; %fixa bias

Interval = quantile(t_th,0.95) % det v?rde som det ?r 95% chans att en v?g ?r lika stort eller mindre ?n

figure
hold on
grid off
cdfplot(t_th);
title('Empirical CDF with Interval at 95% confidence level for a 100-year wave')
line([Interval Interval], [0 1],'Color','red','Linestyle','--');
ylabel('F(100-year wave)');
xlabel('100-year wave (meters)');
hold off
