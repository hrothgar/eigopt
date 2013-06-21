
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example(n)
if nargin < 1, n = 1; end

clear opt;
opt.plotfreq = 25;
opt.dispfreq = 10;
opt.splitting = 1;
opt.endangeredlimit = 1;
opt.cellsmitosisage = Inf;
opt.maxfeval = 2000;
opt.tol = 1e-6;
opt.gamma = -4;

switch n
    case 1
        rng(1);
        pars.A = full(gallery('poisson',15));
        pars.B = randn(225,100);
        [r,bb] = eigopt_multi_mesh('fdist_uncont',[-2 -pi],[6 8-pi],pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 2
        rng(1);
        pars.A = randn(10);
        [r,bb] = eigopt_multi_mesh('fdist_defective',[-2 -pi],[6 8-pi],pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 3
        rng(1);
        pars.A = randn(10);
        [r,bb] = eigopt_multi_mesh('fdist_triple',[-2 -pi],[6 8-pi],pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)

    case 4
        rng(1);
        dim = 3;
        pars.A = randn(10,10,dim+1);
        ee = ones(1,dim);
        [r,bb] = eigopt_multi_mesh('minimize_max_multi',-2*ee,2*ee,pars,opt);
        figure; plot_results(r);
        figure; hist([bb.iternum],15)
end
