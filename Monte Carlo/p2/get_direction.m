function dir = get_direction(r_num, free_neig)
% Finds a new walking direction based on a random number uniform in [0,1]
len = length(free_neig);
index = floor(len*r_num) + 1;
dir = free_neig(index);
end