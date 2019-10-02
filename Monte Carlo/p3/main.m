%% 1c
load('coal_mine_disasters.mat') %T (tau), tidpunkten f?r en olycka
len_tau = length(T);
iter = 1e4;
rho = 1/20;
psi = 1;
d = 3;
burn_in = 2e3;
tau = T;
initial_breakpoints = (1658:322/d:1980);

task1
task1c
for i = 2:d
    mean(t_long(:,i))
end

%%  1d
load('coal_mine_disasters.mat') %T (tau), tidpunkten f?r en olycka
len_tau = length(T);
iter = 1e4;
rho = 1/20;
psi = 0.1;
d = 2;
burn_in = 2e3;
tau = T;
initial_breakpoints = (1658:322/d:1980);
task1
task1d
psi = 1;
task1
task1d
psi = 10;
task1
task1d
%% 1e
task1e

%% task2
task2b
%% 
task2c