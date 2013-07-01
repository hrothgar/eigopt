

opt.plotfreq = Inf;         % plotting frequency (if 2D)
opt.dispfreq = 1;           % verbose output frequency
opt.splitting = 0;          % master flag for splitting
opt.mincellcount = 1;       % split all cells if there are this many or fewer
opt.maxquadspercell = 40;   % maximum number of iterations on a single cell
opt.maxfeval = 4000;        % maximum iteration count, in a sense
opt.tol = 1e-3;             % desired level of accuracy
opt.gamma = -1;             % the gamma thing
opt.minmax = 0;             % minimize (==0) or maximize (~=0) ?

lb = [-2 -4];
ub = [6 4];
func = 'fdist_uncont';

rng(1);
pars.A = full(gallery('poisson',15));
pars.B = randn(225,100);

[r,bb] = eigopt_multi_mesh_breadth(func,lb,ub,pars,opt);
