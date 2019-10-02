%preperation labb1
%% task 1.4
%estimating the params with least squares
%from 3.4
A = [ 1 -1.5 0.7];
C = [1 zeros(1, 11) -0.5];
A5 = [1 zeros(1, 4) -1];
Astar = conv(A, A5);
e = randn(600, 1);
y = filter(C, Astar, e);
y = y(100 : end);
plot(y)


yt = [y(9:end) - y(9-5:end-5)];
%r?tt?
xstar = [e(108:end), e(108-12:end-12), (y(9-1:end-1)-y(9-6:end-6)), (y(9-3:end-3)-y(9-8:end-8))];
%solvs the normal equation
theta = (xstar'*xstar)\(xstar'*yt)

yy = xstar*theta;
figure
plot(yy)
