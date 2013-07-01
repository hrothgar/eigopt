

opts.plotfreq   = 1;      %-- frequency of plotting (if 2-D)
opts.dispfreq   = 1;        %-- frequency of verbose output of status
opts.splitting  = 1;        %-- recursively split into subboxes?
opts.searchtype = 0;        %-- 0 == breadth-first; 1 == depth-first
opts.maxquads   = 40;       %-- maximum number of quadratics per box
opts.maxfeval   = 4000;     %-- maximum number of function evaluations
opts.tol        = 1e-3;     %-- desired level of accuracy
opts.gamma      = -1;       %-- the gamma thing

rng(3);
dim = 4;
clear pars
for d = 1:dim+1,
    A = randn(dim,dim);
    pars.A(:,:,d) = A + A';
end
ee = ones(dim,1);
[fmin, xmin, hist, boxes] = eigopt('minimize_max_multi',[-ee,ee],opts,pars);
