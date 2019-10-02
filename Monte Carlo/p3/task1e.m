%% task1 Gibbs Sampler
load('coal_mine_disasters.mat') %T (tau), tidpunkten f?r en olycka
len_tau = length(T);
iter = 1e4;
d = 3;
tau = T;
psi = 1;

loops = 4; %best?mmer hur m?nga olika rho vi vill kolla 

plotvect_bp_1 = zeros(iter+1,loops);
plotvect_bp_2 = zeros(iter+1,loops);
plotvect_l_1 = zeros(iter+1,loops);
plotvect_l_2 = zeros(iter+1,loops);
plotvect_l_3 = zeros(iter+1,loops);
plotvect_th = zeros(iter+1,loops);

for p = 1:loops
%initial values
theta = gamrnd(2,1/psi);
lambda = gamrnd(2,1/theta, d, 1);
t = (1658:322/d:1980);
n_tau = zeros(d,1);
rho = p/40; 

c = 1;
for i = 1:len_tau
    if tau(i) < t(c+1)
        n_tau(c) = n_tau(c) + 1;
    else
        c = c + 1;
    end
end
t_diff = t(2:end) - t(1:end-1);
th_long = zeros(iter+1,1);
l_long = zeros(iter+1,d);
t_long = zeros(iter+1,d+1);

th_long(1) = theta;
l_long(1,:) = lambda';
t_long(1,:) = t;
% The sick Gibbs Sampler
for k  = 1:iter
    theta = gamrnd(4, 1/sum(lambda) + psi);
    for i = 1:d
        lambda(i) = gamrnd(n_tau(i) + 2, 1/(theta + t_diff(i)));      
    end
    for c = 2:d
        R = rho*(t(c+1)-t(c-1));
        epsilon = (2*R * rand) - R;
        t_star = t(c) + epsilon;
        if t_star > t(c-1) && t_star < t(c+1)
            log_func =@(x, nl, nr)   lambda(c)*(x - t(c+1)) + lambda(c-1)*(t(c-1) - x) ...
            + nr*log(lambda(c)) + nl*log(lambda(c-1)) ...
            + log(t(c+1)*x + x*t(c-1) - x^2 - t(c+1)*t(c-1)) ;
            delta_n = sum(tau<t_star)-sum(tau<t(c));
            nl_star = n_tau(c-1) + delta_n;
            nr_star = n_tau(c) - delta_n;
            alpha = exp(log_func(t_star, nl_star, nr_star) - log_func(t(c), n_tau(c-1), n_tau(c)));
            if alpha > rand
                t(c) = t_star;
                n_tau(c-1) = nl_star;
                n_tau(c) = nr_star;
            end
        end
    end
    t_diff = t(2:end) - t(1:end-1);
    th_long(k+1) = theta;
    l_long(k+1,:) = lambda';
    t_long(k+1,:) = t;

end

%% Mixing

plotvect_bp_1(:,p) = t_long(:,2);
plotvect_bp_2(:,p) = t_long(:,3);
plotvect_l_1(:,p) = l_long(:,1);
plotvect_l_2(:,p) = l_long(:,2);
plotvect_l_3(:,p) = l_long(:,3);
plotvect_th(:,p) = th_long;


end

%% plottar f?r att kolla posteriors
figure
sgtitle('Breakpoint 1')
subplot(411)
hist(plotvect_bp_1(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_bp_1(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_bp_1(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_bp_1(:,4),100);
title('rho = 1/10')

figure
sgtitle('Breakpoint 2')
subplot(411)
hist(plotvect_bp_2(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_bp_2(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_bp_2(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_bp_2(:,4),100);
title('rho = 1/10')

figure
sgtitle('Lambda 1')
subplot(411)
hist(plotvect_l_1(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_l_1(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_l_1(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_l_1(:,4),100);
title('rho = 1/10')

figure
sgtitle('Lambda 2')
subplot(411)
hist(plotvect_l_2(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_l_2(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_l_2(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_l_2(:,4),100);
title('rho = 1/10')

figure
sgtitle('Lambda 3')
subplot(411)
hist(plotvect_l_3(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_l_3(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_l_3(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_l_3(:,4),100);
title('rho = 1/10')

figure
sgtitle('Theta')
subplot(411)
hist(plotvect_th(:,1),100);
title('rho = 1/40')
subplot(412)
hist(plotvect_th(:,2),100);
title('rho = 1/20')
subplot(413)
hist(plotvect_th(:,3),100);
title('rho = 3/40')
subplot(414)
hist(plotvect_th(:,4),100);
title('rho = 1/10')

%% Nytt test f?r att kolla mixingplottar
figure
title('Breakpoint 1')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_bp_1(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_bp_1(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_bp_1(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_bp_1(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off

figure
title('Breakpoint 2')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_bp_2(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_bp_2(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_bp_2(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_bp_2(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off
% 
figure
title('Lambda 1')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_l_1(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_l_1(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_l_1(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_l_1(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off
% 
figure
title('Lambda 2')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_l_2(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_l_2(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_l_2(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_l_2(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off
% 
figure
title('Lambda 3')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_l_3(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_l_3(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_l_3(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_l_3(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off
% 
figure
title('Theta')
xlabel('lag')
ylabel('r(l)')
hold on
plot(autocorr(plotvect_th(:,1),'NumLags',100),'-');
%title('rho = 1/40')
plot(autocorr(plotvect_th(:,2),'NumLags',100),'-');
%title('rho = 1/20')
plot(autocorr(plotvect_th(:,3),'NumLags',100),'-');
%title('rho = 3/40')
plot(autocorr(plotvect_th(:,4),'NumLags',100),'-');
%title('rho = 1/10')
legend('rho = 1/40','rho = 1/20','rho = 3/40','rho = 1/10')
hold off