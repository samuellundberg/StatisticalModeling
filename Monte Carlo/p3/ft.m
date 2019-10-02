function [ Y ] = ft( x, nl, nr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        lz =@(x, nl, nr) log(x) + log(t(c+1) + t(c-1) - x) + (x*(lambda(c)-lambda(c-1))) + nl*log(lambda(c-1)) + nr*log(lambda(c)); 

end

