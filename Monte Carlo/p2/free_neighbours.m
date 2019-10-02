function FN = free_neighbours(m, pos, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dim = length(pos);
free_neig = ones(2*dim,1);
delta = zeros(1, dim);
for i = 1:dim
    delta(i) = 1;
    index = get_index(pos,n,-delta);
    if m(index)
       free_neig(2*i - 1) = 0;
    end
    index = get_index(pos,n,delta);
    if m(index)
       free_neig(2*i) = 0;
    end
    delta(i) = 0;
end
FN = find(free_neig);

end
