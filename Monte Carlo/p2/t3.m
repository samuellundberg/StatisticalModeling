
N = 1e3;    % particles
n = 20;     % steps
dim = 2;      % dimensions
success = ones(N, n);
gn = ones(N, n);
walks = rand(N,n);
pos = (n + 2)*ones(N,dim);
meshes = cell(N,1);
nbr_of_elements = (2*(n+2) -1 )^dim;
temp_v = sparse(nbr_of_elements,1);
index = get_index(pos(1,:),n);
temp_v(index) = 1;
for k = 1:N
    meshes{k} = temp_v;
end
for i = 1:n
    for k = 1:N
        if success(k, i) == 0
            continue
        end
        gn(k,i:n) = gn(k,i:n)/4;
        dir = floor(4*walks(k, i)) + 1;
        pos(k,:) = update_pos(pos(k,:), dir);
        index = get_index(pos(k,:), n);
        if meshes{k}(index)
            success(k, i:end) = 0;
        else 
           meshes{k}(index) = 1;
        end
    end
end
% makes w as a weight function for a successfull walk
g_mod = ones(N,n);
for i = 1:n
    g_mod(:,i:n) = g_mod(:,i:n)/4;
end
success_ratio = sum(success) ./ N;
%c_simp = success_ratio ./ g
c_simp = sum(success ./ (N*g_mod))';
individual_esimates = success ./ (g_mod);
interval = (1.96 * std(individual_esimates) /sqrt(N))';
conf_int = [c_simp - interval, c_simp + interval];

[(1:5)', conf_int(1:5,:)]
[10, conf_int(10,:)] 
[20, conf_int(20,:)]




