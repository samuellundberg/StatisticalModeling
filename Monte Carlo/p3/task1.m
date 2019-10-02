%% task1 Gibbs Sampler
%initaial values
theta = gamrnd(2,1/psi);
lambda = gamrnd(2,1/theta, d, 1);
t = initial_breakpoints;
n_tau = zeros(d,1);
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
    theta = gamrnd(4, 1/(sum(lambda) + psi));
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
th_long = th_long(burn_in+1:end);
l_long = l_long(burn_in+1:end, :);
t_long = t_long(burn_in+1:end, :);