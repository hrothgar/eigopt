function [f, g] = minimize_max_multi(omega,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fdist_instab(omega,pars)
% task:
%		calculates f = sigma_min(A - omega*i*I) and its
%		derivative g
% note:
%		the input matrix A must be passed through pars.A

MA = pars.A;
[m,n,p] = size(MA);
l = length(omega);

if ((n ~= m) | (p-1 ~= l))
	error('the input matrix must be square');
end


% keyboard
A = MA(:,:,1);
for j = 2:p
	A = A + omega(j-1,1)*MA(:,:,j);
end



if (n <= 150)


	[V,D] = eig(A);

%	diag(D)
	f = D(n,n);


	for j = 1:p-1
		g(j,1) = real(V(:,n)'*MA(:,:,j+1)*V(:,n));
	end
else
			 
	[eigvec,d] = eigs(A,1,'lm');

	vm = eigvec(1:n,1);
	vm = vm/norm(vm);

	f = real(d(1,1));
	
	for j = 1:p-1
		g(j,1) = real(V(:,n)'*MA(:,:,j)*V(:,n));
	end				  

end




return;