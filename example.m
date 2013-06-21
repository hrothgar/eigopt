
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example(n)
if nargin < 1, n = 1; end


%%% OUTPUT
opt.plotfreq = 20;          % plotting frequency (if 2D)
opt.dispfreq = 1;           % verbose output frequency

%%% MESH
% opt.splitting = 1;          % master flag for splitting
% opt.depthfirst = 0;         % flag for search type; 0 => breadth-first
opt.mincellcount = 1;       % split all cells if there are this many or fewer
opt.maxquadspercell = 80;   % maximum number of iterations on a single cell

%%% STOPPING CRITERIA
opt.maxfeval = 2000;        % maximum iteration count, in a sense
opt.tol = 1e-6;             % desired level of accuracy

%%% THE PROBLEM
opt.gamma = -4;             % the gamma thing
opt.minmax = 0;             % minimize (==0) or maximize (~=0) ?


switch n
    case 1
        rng(1);
        pars.A = full(gallery('poisson',15));
        pars.B = randn(225,100);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_depth('fdist_uncont',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 2
        rng(1);
        pars.A = randn(10);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_depth('fdist_defective',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 3
        rng(1);
        pars.A = randn(10);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_depth('fdist_triple',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 4
        rng(1);
        dim = 5;
        for d = 1:dim+1,
            A = randn(dim,dim);
            pars.A(:,:,d) = A + A';
        end
        ee = ones(1,dim);
        [r,bb] = eigopt_multi_mesh_depth('minimize_max_multi',-2*ee,2*ee,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)
end
