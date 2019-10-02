function [pos] = update_pos(pos, dir)

a = mod(dir,2);
if ~a
    a = -1;
end
b = round(dir/2);
pos(b) = pos(b) - a;

end
