
rng(1);
dim = 3;
clear expars;
for d = 1:dim+1,
    A = randn(dim,dim)/2;
    expars.A(:,:,d) = A + A';
end
bounds = [-ones(dim,1) ones(dim,1)];
opts.plotfreq = Inf;
opts.tol = 1e-3;
opts.maxquads = 40;
opts.dispfreq = 1;
opts.gamma = -0.1;

[fmin,xmin,hist,boxes] = eigopt('minimize_max_multi', bounds, opts, expars);

% figure; plot_results(r);
% figure; hist([bb.iternum],15);
% title('Histogram of the age of mesh cells')
% xlabel('iterations')