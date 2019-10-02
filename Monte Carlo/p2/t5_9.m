%% t5-9. 5&6 id d = 2, otherwise it is 9
N = 1e4;    % particles
n = 20;     % steps
dim = 2;      % dimensions
success = ones(N, n);
gn = ones(N, 1);
walks = rand(N,n);
pos = (n + 2)*ones(N,dim);
meshes = cell(N,1);
nfn = zeros(N,n);
nbr_of_elements = (2*(n+2) -1 )^dim;
temp_v = sparse(nbr_of_elements,1);
index = get_index(pos(1,:),n);
temp_v(index) = 1;

c_n = ones(n,1);
sigma_cn = ones(n,1);

for k = 1:N
    meshes{k} = temp_v;
end
 % loop over steps
for i = 1:n
    % loop over particles
    for k = 1:N
        free_neig = free_neighbours(meshes{k}, pos(k,:), n);
        if isempty(free_neig)
            success(k, i) = 0;
            gn(k) = inf;
            continue
        end
        nfn(k,i) = length(free_neig);
        dir = get_direction(walks(k,i), free_neig);
        gn(k) = 1 / length(free_neig);
        pos(k,:) = update_pos(pos(k,:), dir);
        index = get_index(pos(k,:),n);
        meshes{k}(index) = 1;
        
    end
    
     if i == 1
        c_n(i) = sum(1./ (N*gn));
        sigma_cn(i) = std(1./ (1*gn));
     else
        c_n(i) = c_n(i-1)*sum(1./ (N*gn));
        sigma_cn(i) = std(c_n(i-1)./ (1*gn));
     end
     new_particles = randsample(N,N, true, 1./ (N*gn));
     for k = 1:N
         new_p = new_particles(k);
         pos(k,:) = pos(new_p,:);
         meshes{k} = meshes{new_p};
         % ska nfn updateras h?r
     end
end

success_ratio = sum(success) ./ N;
c_n;


interval = (1.96 * sigma_cn /sqrt(N));
conf_int = [c_n - interval, c_n + interval];
if n > 4
    [(1:5)', conf_int(1:5,:)]
end
if n > 19
    [10, conf_int(10,:)] 
    [20, conf_int(20,:)]
end
%% t6

cut_off = 0;
N, n, dim, cut_off
gamma = zeros(3,1);
mu = zeros(3,1);
A = zeros(3,1);
% Endast A som ok?nnd
gamma(1) = 43/32;
mu(1) = mean(nfn(:));
X = ones(n-cut_off,1);
y = log(c_n(cut_off+1:end))- log((cut_off+1:n))'.*(gamma(2)-1) - (cut_off+1:n)'.*log(mu(1));
A(1) = exp(X\y)

% A och mu ok?nnda NOT USED
gamma(2) = 43/32;
X = [ones(n-cut_off,1), (cut_off+1:n)'];
y = log(c_n(cut_off+1:end)) - log((cut_off+1:n))'.*(gamma(2)-1);
omega = X\y;
A(2) = exp(omega(1));
mu(2) = exp(omega(2));

% alla ?r ok?nda
X = [ones(n-cut_off,1), (cut_off+1:n)', log((cut_off+1:n)')];
y = log(c_n(cut_off+1:end));
omega = X\y;
A(3) = exp(omega(1));
mu(3) = exp(omega(2));
gamma(3) = omega(3) + 1;

[gamma mu A]

