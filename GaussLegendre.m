function [x, w] = GaussLegendre(n)
i   = 1:n-1;
a   = i./sqrt(4*i.^2-1);
CM  = diag(a,1) + diag(a,-1);
[V L]   = eig(CM);
[x ind] = sort(diag(L));
V       = V(:,ind)';
w       = 2 * V(:,1).^2;
