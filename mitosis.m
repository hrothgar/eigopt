% divides a box into 2^dim smaller cells of equal size, e.g.
% > b0 = [0 0];
% > b1 = [1 1];
% > [x0, x1] = mitosis(b0, b1)
% x0, x1 have 2^dim rows, each corresponding to the a new cell

function [x0, x1] = mitosis(b0, b1)

dim = length(b0);
weights = de2bi((0:2^dim-1)')/2;
dx = ones(2^dim,1)*(b1(:)-b0(:))';

x0 = ones(2^dim,1)*b0(:)' + weights.*dx;
x1 = x0 + dx/2;

return
