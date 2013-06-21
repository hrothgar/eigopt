function [f, g] = minimize_quad_max_multi(omega,pars)


A0 = pars.A0;
A1 = pars.A1;
A2 = pars.A2;
A3 = pars.A3;

[m,n] = size(A0);


[V,D] = eig(A0+omega(1)^2*A1+omega(2)^2*A2+omega(1)*omega(2)*A3);

f = D(n,n);
g(1) = V(:,n)'*(2*omega(1)*A1 + omega(2)*A3)*V(:,n);
g(2) = V(:,n)'*(2*omega(2)*A2 + omega(1)*A3)*V(:,n);



return;