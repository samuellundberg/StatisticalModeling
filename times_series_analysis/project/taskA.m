%task a. extrapolate the rain with help from kalman
%MODEL THE RAIN AS A AR(1) PROCESS AND REC. USING KALMAN
load proj18
% THE SUM OF THE RAIN SHOULD BE THE SAME FOR OUR METHOD AS THE GIVEN ONE
rain = ElGeneina.rain_org;
rain = sqrt(rain);
MEAN_VALUE_RAIN = mean(rain)
rain = rain - MEAN_VALUE_RAIN;
plot(rain)
hold on
my_plotter(rain)
%% RECONSTRUCT THE RAIN xt FROM THE RAIN yt, WE HAVE TO ASSUME a I KNOWN
y = rain;
N = length(y);
A = 1;
Re = 10; %1
Rw = 20; %5
C = 3; %not time dependent
Rxx1 = Rw;
ztt1 = 0; %z = [xt]
%Vector to store values in
zsave = zeros(1, 3*N);
% Kalman filter. We loop througt y and get three x each 
% time through a nested loop
for k = 1:N
    for kk = 1:3
        Ryy = C*Rxx1*C' + Rw; 
        Kt = Rxx1*C'/Ryy;
        ztt = ztt1 + Kt*(y(k) - C*ztt1);
        Rxx = (eye(length(C)) - Kt * C)* Rxx1;
        %Save
        zsave(:, 3*(k-1) + kk) = ztt;
        %Predict
        Rxx1 = A * Rxx*A' + Re;
        ztt1 = A*ztt;
    end
end
figure
plot(zsave, '*-')
%verifying that the sum of the rain is right 
scalefactor = sum(y + MEAN_VALUE_RAIN)/sum(zsave(15:end) + MEAN_VALUE_RAIN/3)
%%

zsum = zeros(size(y));
for i=0:length(y)-3
    zsum(i+1) = sum(zsave(3*i + 1:3*i + 3));
end
figure
subplot(2,1,1)
plot(zsum)
title('Original rain compared to the truncated reconstruction')
hold on
plot(y, 'r-.')
legend('Truncated reconstruction', 'Original rain')

subplot(2,1,2)
plot([39:85],zsum(39:85))
title('A closer look of the same thing')

hold on
plot([39:85],y(39:85), 'r-.')
axis tight

clear y N a1 A Re Rw C Rxx1 ztt1 Ryy Kt ztt Rxx
clear scalefactor zsum kk k i
