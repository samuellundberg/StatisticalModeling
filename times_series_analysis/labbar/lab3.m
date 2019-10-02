load tar2.dat
load thx.dat
figure
subplot(2,1,1)
plot(tar2);
subplot(2,1,2)
plot(thx);      %parametrar i tar2. a1 ?ndras
hold on
%% q1
%what is nb and nk b_-pol and kal
na = 2;
model = [na];%,nb,nk]
lambda = 0.9424; %0.95 ives a good estimate with reasonable variance
%use recursiveARX instead of rarx on newer matlab
[Aest, yhat, covAest, yprev] = rarx(tar2, model, 'ff', lambda);

% figure
subplot(2,1,1)
plot(yhat);
subplot(2,1,2)
plot(Aest(:,1), 'b');
hold on
plot(Aest(:,2), 'r');
%% q2
n = 100;
lambdaline = linspace(0.85, 1, n);
ls2 = zeros(n, 1);
for i = 1:length(lambdaline)
[Aest, yhat, CovAest, trash] = rarx(tar2, [2], 'ff', lambdaline(i));
ls2(i) = sum((tar2 - yhat).^2);
end
plot(lambdaline, ls2) 
hold on
argmin = find((ls2-min(ls2))<=0)
opt_lambda = lambdaline(argmin)
plot(opt_lambda, ls2(argmin),'*')

%% q3 we should estimate the parameters instead. see page 300
%Example of Kalman filter
%Simulate process
y = tar2;
% Data length
N = length(y);
% Define the state space equations
A = eye(2); %[-thx(1,:); 1 0]';       %initial
%e_t, model imperfection
Re = [1 0; 0 0]/100; %Hidden state noise covariance matrix
Rw = 125; %Observation variance, m?tfelet, angiver i uppgiften
%usually C should be set here to,
%but in this case C is a function of time .
%Set some initial values
Rxx1 = Rw * eye(2); %Initial variance
xtt1 = [thx(1,:)]'; %Initial state (m0)
%Vector to store values in
xsave = zeros(2, N);
%xsave(:, 1) = Aest(1, :);
% Kalman filter. Start from k=3,
%since we need old values of y.
for k = 3:N
    %Cis, in our case, a function of time.
    C = [-y(k-1), -y(k-2)]; 
    % Update
    Ryy = C*Rxx1*C' + Rw; %fr?n 8.10
    Kt = Rxx1*C'/Ryy;  %fr?n 8.10
    xtt = xtt1 + Kt*(y(k) - C*xtt1);%fr?n 8.10
    Rxx = (eye(2) - Kt * C)* Rxx1;  %fr?n 8.10
    %Save
    xsave(:, k) = xtt;
    %Predict
    Rxx1 = A * Rxx*A' + Re;
    xtt1 = A*xtt;
end;
figure
plot(xsave(1,:), 'b')
hold on
plot(xsave(2,:),'r')

%% q4 vad ?r b?st
rls_err = thx - Aest;
kal_err = thx - xsave';
rls_err = rls_err(100:end,:);
kal_err = kal_err(100:end,:);

whitenessTest(rls_err(:,2))
my_plotter(rls_err(:,2))
figure
whitenessTest(kal_err(:,2))
my_plotter(kal_err(:,2))

%% q5

R = randi(8,500,1);
u = zeros(length(R), 1)';
for i = 2:length(R)
    if u(i-1) && R(i) > 1
        u(i) = 1;
    elseif ~u(i-1) && R(i) == 1
        u(i) = 1;
    end
end
xfake = filter(1, [1 -1], randn(1,500));
b = 20;
Rw = 4;
y = xfake + b*u + Rw*randn(1,500);
plot(y)
%%
% Data length
N = length(u);
% Define the state space equations
A = [1 0; 0 1]; %vi vill inte ?ndra b

%e_t, model imperfection
Re = 0.01*[1 0; 0 0]; %Hidden state noise covariance matrix
Rw = 1; %Observation variance, m?tfelet, angiver i uppgiften
%usually C should be set here to,
%but in this case C is a function of time .
%Set some initial values
Rxx1 = Rw*eye(2);% * eye(2); %Initial variance
ztt1 = [0 b]'; %z = [x, b]'
%Vector to store values in
xsave = zeros(2, N);
%xsave(:, 1) = Aest(1, :);
% Kalman filter. Start from k=3,
%since we need old values of y.
for k = 1:N
    %Cis, in our case, a function of time.
    C = [1, u(k)]; 
    Ryy = C*Rxx1*C' + Rw; %fr?n 8.10
    Kt = Rxx1*C'/Ryy;  %fr?n 8.10
    ztt = ztt1 + Kt*(y(k) - C*ztt1);%fr?n 8.10
    Rxx = (eye(2) - Kt * C)* Rxx1;  %fr?n 8.10
    %Save
    zsave(:, k) = ztt;
    %Predict
    Rxx1 = A * Rxx*A' + Re;
    ztt1 = A*ztt;
end

figure
plot(zsave(1,:))
hold on
plot(xfake, '-.r')
figure
plot(zsave(2,:))
hold on
plot([0 500],[20 20], '-.r')








