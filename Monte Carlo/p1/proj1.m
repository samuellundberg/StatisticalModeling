%% TASK 2
load('powercurve_V90.mat')
points = 100000;
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2 2 2 1.9 1.9 1.9 1.9 1.9 2 1.9 2 2];
weibull = zeros(points,12);
for i = 1:points
    weibull(i,:) = wblrnd(lambda, k);
end
power_curve = zeros(12,1);
conf_factor = zeros(12,1);
for i = 1:12
    v = weibull(:,i);
    power_curve(i) = mean(P(v));
    conf_factor(i) = 1.96*std(P(v))/sqrt(points);
end
conf_int = [power_curve - conf_factor, power_curve + conf_factor]
%% del tva
v_t = zeros(12,points);
pc_t = zeros(12,1);
cf_t = zeros(12,1);
u = rand(points,1);
for i = 1:12
    F4  = wblcdf(4, lambda(i), k(i));
    F25 = wblcdf(25, lambda(i), k(i));
    ustar = u*(F25-F4) + F4;
    x = wblinv(ustar, lambda(i), k(i));
    v_t(i,:) = x;
    pc_t(i) = mean(P(v_t(i,:)))*(F25-F4);
    cf_t(i) = 1.96*std(P(v_t(i,:)))/sqrt(points);
end
ci_t = [pc_t - cf_t, pc_t + cf_t]
%% 2b importance sampling, g kan f?rb?ttras
pc_is = zeros(12,1);
cf_is = zeros(12,1);
u_1     = u;
sigma = 4.5;
mu = 1.25;
for i = 1:12
    x_is  = norminv(u_1, mu*lambda(i), sigma);
    f     = wblpdf(x_is, lambda(i), k(i));
    g     = normpdf(x_is, mu*lambda(i), sigma);
    phi   = P(x_is);
    omega = f./g;
    PO    = phi.*omega;
    pc_is(i)   = mean(PO);
    cf_is(i)   = 1.96 * std(PO) / sqrt(points);
end
ci_is = [pc_is - cf_is, pc_is + cf_is]
%% 2c Antithetic sampling
pc_as = zeros(12,1);
cf_as = zeros(12,1);
u_temp     = u(1:end/2);
W = zeros(points/2, 12);
for i = 1:12
    F4  = wblcdf(4, lambda(i), k(i));
    F25 = wblcdf(25, lambda(i), k(i));
    u_star = u_temp*(F25-F4) + F4;
    usp = (1-u_temp)*(F25-F4) + F4;
    x = wblinv(u_star, lambda(i), k(i));
    x_prim = wblinv(usp, lambda(i), k(i));
    V = P(x);
    V_prim = P(x_prim);
    W(:,i) = (V+V_prim)/2;
    pc_as(i) = mean(W(:,i))*(F25-F4);
    cf_as(i) = 1.96*std(W(:,i))/sqrt(points/2);
end
ci_as = [pc_as - cf_as, pc_as + cf_as]
%% 2d Estimate , p(P(V) > 0).
% P(V) > 0|4<=v<=25 = P(25)-P(4)
p = wblcdf(25, lambda, k) - wblcdf(4, lambda, k);
mean(p)

%% 2e calc E(P)/E(Ptot) m.h.a. bra wind use CMC???
rho = 1.225;
d = 90;
mean_Ptot = rho * pi * d^2 * (lambda.^3 .* gamma(1 + 3./k)) / 8;
ci_wr= zeros(12,2);
for i = 1:12
    ci_wr(i,:) = ci_as(i,:)/mean_Ptot(i);
end
ci_wr
%% 2f
capacity_factor     = pc_as/3000000;        % MEDELV?RDE
Availability_factor = p;                    % fr?n 2d
mean(capacity_factor)
mean(Availability_factor)

%% TASK 3

load('powercurve_V90.mat')
points = 10000;
lambda = 9.13;
rho = 1.225;
d = 90;
k = 1.96;
a = 0.638;
p = 3;
q = 1.5;

F = @(v_1,v_2) wblcdf(v_1,lambda,k).*wblcdf(v_2,lambda,k).*(1+a.*(1-wblcdf(v_1,lambda,k).^p).^q.*(1-wblcdf(v_2,lambda,k).^p).^q);
f = @(v_1,v_2) wblpdf(v_1,lambda,k).*wblpdf(v_2,lambda,k).*(1+a.*(1-wblcdf(v_1,lambda,k).^(p)).^(q-1).*(1-wblcdf(v_2,lambda,k).^(p)).^(q-1).*(wblcdf(v_1,lambda,k).^(p).*(1+p.*q)-1).*(wblcdf(v_2,lambda,k).^(p).*(1+p.*q)-1));
omega_fkt = @(v_1,v_2) (1+a.*(1-wblcdf(v_1,lambda,k).^(p)).^(q-1).*(1-wblcdf(v_2,lambda,k).^(p)).^(q-1).*(wblcdf(v_1,lambda,k).^(p).*(1+p.*q)-1).*(wblcdf(v_2,lambda,k).^(p).*(1+p.*q)-1));

g_new = @(v_1,v_2) wblpdf(v_1,lambda,k).*wblpdf(v_2,lambda,k);
P_new = @(v_1,v_2) P(v_1)+P(v_2);

g_2_new = zeros(points,1);
phi_new_1 = zeros(points,1);
phi_new_2 = zeros(points,1);
f_2_new = zeros(points,1);
omega_new = zeros(points,1);

P_1 = zeros(points,1);
P_2 = zeros(points,1);

wind_1 = wblrnd(lambda,k,points,1);
wind_2 = wblrnd(lambda,k,points,1);
%%
for i = 1:points
        P_1(i) = P(wind_1(i));
        P_2(i) = P(wind_2(i));
    
        f_2_new(i) = f(wind_1(i),wind_2(i));  
        g_2_new(i) = g_new(wind_1(i),wind_2(i));
        phi_new_1(i) = P(wind_1(i))+P(wind_2(i)); %f?r a-delen
        phi_new_2(i) = P(wind_1(i)).*P(wind_2(i)); %f?r b-delen
        omega_new(i) = omega_fkt(wind_1(i),wind_2(i));     
end
        PO_new_1 = phi_new_1.*omega_new; %a-del
        PO_new_2 = phi_new_2.*omega_new; %b-del
        P_1_new = P_1.*omega_new;
        P_2_new = P_2.*omega_new;
%% 3a 1d because they are independent
expvalue = mean(PO_new_1)
%% 3b
Cov_term_1 = mean(PO_new_2); 
Cov_term_2 = mean(P_1_new)*mean(P_2_new);
covariance = Cov_term_1 - Cov_term_2
%% 3c
variab = var(P_1_new) + var(P_2_new) + 2*covariance
st_dev = sqrt(variab)
%% 3d we use importance sampling 
phi_new_1_d = zeros(points,1);
phi_new_2_d = zeros(points,1);
for i = 1:points
        phi_new_1_d(i) = P(wind_1(i))+P(wind_2(i)) > 3*10^6; %st?rre ?n
        phi_new_2_d(i) = P(wind_1(i))+P(wind_2(i)) < 3*10^6; %mindre ?n   
end
PO_new_1_d = phi_new_1_d.*omega_new; %st?rre ?n
PO_new_2_d = phi_new_2_d.*omega_new; %mindre ?n 
mean(PO_new_1_d)
interval_1 = mean(PO_new_1_d) + (1.96*(std(PO_new_1_d)))/sqrt(points)*[-1 1]
mean(PO_new_2_d)
interval_2 = mean(PO_new_2_d) + (1.96*(std(PO_new_2_d)))/sqrt(points)*[-1 1]
% does not quite sum to one