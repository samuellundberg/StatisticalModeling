N = 1e4;    % particles
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
% loop over steps
for i = 1:n
    % loop over particles
    for k = 1:N
        if success(k, i) == 0
            continue
        end
        free_neig = free_neighbours(meshes{k}, pos(k,:), n);
        if isempty(free_neig)
            success(k, i:end) = 0;
            gn(k, i:n) = inf;
            continue
        end
        dir = get_direction(walks(k,i), free_neig);
        gn(k, i:n) = gn(k, i:n) / length(free_neig);
        pos(k,:) = update_pos(pos(k,:), dir);
        index = get_index(pos(k,:),n);
        meshes{k}(index) = 1;
    end
end
success_ratio = sum(success) ./ N;
% c = ratio / g
c_simp = sum(1./(N*gn))';

individual_esimates = 1 ./ (gn);
interval = (1.96 * std(individual_esimates) /sqrt(N))';
conf_int = [c_simp - interval, c_simp + interval];

[(1:5)', conf_int(1:5,:)]
[10, conf_int(10,:)] 
[20, conf_int(20,:)]
