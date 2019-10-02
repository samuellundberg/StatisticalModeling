%task 2
load('powercurve_V90.mat')
points = 10000;
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2 2 2 1.9 1.9 1.9 1.9 1.9 2 1.9 2 2];
weibull = zeros(points,12);
for i = 1:points
    weibull(i,:) = wblrnd(lambda, k);
end
power_curve = zeros(12,1);
conf_factor = zeros(12,1);
v_t = zeros(12,points);
pc_t = zeros(12,1);
cf_t = zeros(12,1);
pc_is = zeros(12,1);
cf_is = zeros(12,1);
pc_as = zeros(12,1);
cf_as = zeros(12,1);

conf_int = zeros(100,1);
ci_t  = zeros(100,1);
ci_is = zeros(100,1);
ci_as = zeros(100,1);

u = rand(points,1);
u_temp     = u(1:end/2);
sigma = 4.5;
mu = 1.25;
W = zeros(points/2, 12);


for i = 1:12
    v = weibull(:,i);
    power_curve(i) = mean(P(v));
    conf_factor(i) = 1.96*std(P(v))/sqrt(points);
    
    F4  = wblcdf(4, lambda(i), k(i));
    F25 = wblcdf(25, lambda(i), k(i));
    ustar = u*(F25-F4) + F4;
    x = wblinv(ustar, lambda(i), k(i));
    v_t(i,:) = x;
    pc_t(i) = mean(P(v_t(i,:)))*(F25-F4);
    cf_t(i) = 1.96*std(P(v_t(i,:)))/sqrt(points);
    
    x_is  = norminv(u, mu*lambda(i), sigma);
    f     = wblpdf(x_is, lambda(i), k(i));
    g     = normpdf(x_is, mu*lambda(i), sigma);
    phi   = P(x_is);
    omega = f./g;
    PO    = phi.*omega;
    pc_is(i)   = mean(PO);
    cf_is(i)   = 1.96 * std(PO) / sqrt(points);
    
    u_star = ustar(1:end/2);
    usp = (1-u_temp)*(F25-F4) + F4;
    x = wblinv(u_star, lambda(i), k(i));
    x_prim = wblinv(usp, lambda(i), k(i));
    V = P(x);
    V_prim = P(x_prim);
    W(:,i) = (V+V_prim)/2;
    pc_as(i) = mean(W(:,i))*(F25-F4);
    cf_as(i) = 1.96*std(W(:,i))/sqrt(points/2);
    
end
conf_int = [mean(power_curve) - mean(conf_factor), mean(power_curve) + mean(conf_factor)]
ci_t  = [mean(pc_t) - mean(cf_t), mean(pc_t) + mean(cf_t)]
ci_is = [mean(pc_is) - mean(cf_is), mean(pc_is) + mean(cf_is)]
ci_as = [mean(pc_as) - mean(cf_as), mean(pc_as) + mean(cf_as)]




plot(1,conf_int(1), 'm')
hold on
plot(1,conf_int(2), 'm')
plot(1, ci_t(1), 'g')
plot(1, ci_t(2), 'g')
plot(1, ci_is(1), 'r')
plot(1, ci_is(2), 'r')
plot(1, ci_as(1), 'b')
plot(1, ci_as(2), 'b')
