
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example_breadth(n, dimension)
if nargin < 1, n = 1; end
if nargin < 2, dimension = 2; end


%%% OUTPUT
opt.plotfreq = 20;          % plotting frequency (if 2D)
opt.dispfreq = 1;           % verbose output frequency

%%% MESH
opt.splitting = 1;          % master flag for splitting
opt.mincellcount = 1;       % split all cells if there are this many or fewer
opt.maxquadspercell = 80;   % maximum number of iterations on a single cell

%%% STOPPING CRITERIA
opt.maxfeval = 2000;        % maximum iteration count, in a sense
opt.tol = 1e-3;             % desired level of accuracy

%%% THE PROBLEM
opt.gamma = -1;             % the gamma thing
opt.minmax = 0;             % minimize (==0) or maximize (~=0) ?


switch n
    case 1
        rng(1);
        pars.A = full(gallery('poisson',15));
        pars.B = randn(225,100);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_breadth('fdist_uncont',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15);
        title('Histogram of the age of mesh cells')
        xlabel('iterations')

    case 2
        rng(1);
        pars.A = randn(10);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_breadth('fdist_defective',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15);
        title('Histogram of the age of mesh cells')
        xlabel('iterations')

    case 3
        rng(1);
        pars.A = randn(10);
        lb = [-2 -pi];
        ub = [6 8-pi];
        [r,bb] = eigopt_multi_mesh_breadth('fdist_triple',lb,ub,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15);
        title('Histogram of the age of mesh cells')
        xlabel('iterations')

    case 4
        rng(1);
        dim = dimension;
        for d = 1:dim+1,
            A = randn(dim,dim);
            pars.A(:,:,d) = A + A';
        end
        ee = ones(1,dim);
        [r,bb] = eigopt_multi_mesh_breadth('minimize_max_multi',-ee,ee,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15);
        title('Histogram of the age of mesh cells')
        xlabel('iterations')
end
