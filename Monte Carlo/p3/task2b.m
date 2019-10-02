%% b
load('atlanticdata.mat')
F_inv = @(u,mu,beta) mu-beta.*(log(-log(u)));
atlantic = X;
l = length(atlantic);
n = 20000;
[beta mu] = est_gumbel(atlantic);

gumbel_est = zeros(n,2);
t = zeros(n,2);
delta_Y = zeros(n,2);

for i = 1:n
    P0_hat = F_inv(rand(l,1),mu,beta);
    [x y] = est_gumbel(P0_hat);   
    gumbel_est(i,1) = x;
    gumbel_est(i,2) = y;
    t(i,1) = gumbel_est(i,1)-beta; %detta g?r vi f?r att fixa bias
    t(i,2) = gumbel_est(i,2)-mu;
end


delta_beta_mean = mean(t(:,1)); % detta g?r vi f?r att fixa bias
delta_mu_mean = mean(t(:,2));

gumbel_est(:,1) = gumbel_est(:,1)-delta_beta_mean;
gumbel_est(:,2) = gumbel_est(:,2)-delta_mu_mean;

I_mu = [quantile(gumbel_est(:,2),0.025) quantile(gumbel_est(:,2),0.975)]
I_beta = [quantile(gumbel_est(:,1),0.025) quantile(gumbel_est(:,1),0.975)]

figure
hold on
grid off
subplot(2,1,1);
cdfplot(gumbel_est(:,1));
title('Empirical CDF with Interval at 95% confidence level for beta')
line([I_beta(1) I_beta(1)], [0 1],'Color','red','Linestyle','--');
line([I_beta(2) I_beta(2)], [0 1],'Color','red','Linestyle','--');
ylabel('F(Beta)');
xlabel('Beta');

subplot(2,1,2);
cdfplot(gumbel_est(:,2));
title('Empirical CDF with Interval at 95% confidence level for mu')
line([I_mu(1) I_mu(1)], [0 1],'Color','red','Linestyle','--');
line([I_mu(2) I_mu(2)], [0 1],'Color','red','Linestyle','--');
ylabel('F(mu)');
xlabel('mu');
hold off

