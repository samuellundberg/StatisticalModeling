% Simulate ut in 3.3 with the markov chain
% Estimate the parameters of an AR(2)-process using a Kalman filter. 
% Section 4 gives a rough outline, which might be of some help.
% yt = xt + b*ut + vt <=> ut = (yt - xt - vt)/b
R = randi(8,500,1);
u = zeros(length(R), 1);
for i = 2:length(R)
    if u(i-1) && R(i) > 1
        u(i) = 1;
    elseif ~u(i-1) && R(i) == 1
        u(i) = 1;
    end
end
plot(u)

%%
%Example of Kalman filter
%Simulate process
y = u;
% Data length
N = length(y);
% Define the state space equations
A = [7 1; 1 7]/8;       %markov
Re = [1 1; 1 1]; %Hidden state noise covariance matrix
Rw = 1; %Observation variance
%usually C should be set here to,
%but in this case C is a function of time .
%Set some initial values
Rxx1 = 1 * eye(2); %Initial variance
xtt1 = [1 1]'; %Initial state
%Vector to store values in
xsave = zeros(2, N);
% Kalman filter. Start from k=3,
%since we need old values of y.
for k = 3:N
    %Cis, in our case, a function of time.
    C = [k k];
    % Update
    Ryy = 1;
    Kt = 1;
    xtt = 1;
    Rxx = 1;
    %Save
    xsave(:, k) = [1 k];
    %Predict
    Rxx1 = 1;
    xtt1 = 1;
end;

subplot(2,1,1)
plot(xsave(1,:))
subplot(2,1,2)
plot(xsave(2,:))
