function [ index ] = get_index(pos, n, delta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    dim = length(pos);
    p_vec = pos - [0, ones(1,dim-1)];
    if nargin == 3
        p_vec = p_vec + delta;
    end    
    s_vec = (2*(n+2)-1).^(0:dim-1)';
    index = p_vec*s_vec;
end
