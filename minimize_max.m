function [f, g] = minimize_max(omega,pars)




A0 = pars.A0;
A1 = pars.A1;
[n,m] = size(A0);
[p,q] = size(A1);


if (n ~= m)
	error('the input matrix must be square');
end

if ((n ~= p) | (m ~= q))
	error('A0 and A1 must be of same size');
end



if (n <= 150)
	[V,D] = eig(A0 + omega*A1);

	f = D(n,n);
	g = real(V(:,n)'*A1*V(:,n));
else
	[eigvec,d] = eigs(A0 + omega*A1,1,'lm');

	vm = eigvec(1:n,1);
	vm = vm/norm(vm);

	f = real(d(1,1));
	g = real(vm'*A1*vm);
			 

end




return;