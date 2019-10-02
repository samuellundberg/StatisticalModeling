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