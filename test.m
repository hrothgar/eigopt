
rng(1);
pars.A = full(gallery('poisson',15));
pars.B = randn(225,100);

funname = 'fdist_uncont';
lb = [-2 -pi];
ub = [6 8-pi];

clear opt;
opt.plotfreq = Inf;
opt.dispfreq = Inf;
opt.splitting = 1;
opt.endangeredlimit = 1;
opt.cellsmitosisage = Inf;
opt.maxfeval = 2000;
opt.tol = 1e-6;
opt.gamma = -4;

tic;
[r,bb] = eigopt_multi_mesh(funname,lb,ub,pars,opt);
t1 = toc;

tic;
[f, z, lbound, nfevals] = eigopt_multi_main(funname, lb, ub, ...
                    opt.gamma, opt.maxfeval, opt.tol, pars, [], 0);
t2 = toc;

disp(['Mesh: f=' num2str(r(end).f) ', #evals=' num2str(r(end).nfevals) ', t=' num2str(t1)]);
disp(['Main: f=' num2str(f) ', #evals=' num2str(nfevals) ', t=' num2str(t2)]);
