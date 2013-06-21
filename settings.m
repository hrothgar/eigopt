
clear opt;

%%% OUTPUT
opt.plotfreq = 20;          % plotting frequency (if 2D)
opt.dispfreq = 1;           % verbose output frequency

%%% MESH
opt.splitting = 1;          % master flag for splitting
opt.depthfirst = 1;         % flag for search type; 0 => breadth-first
opt.mincellcount = 1;       % split all cells if there are this many or fewer
% opt.maxcellcount = 100;      % maximum number of cells at any point
opt.maxquadspercell = 40;   % maximum number of iterations on a single cell

%%% STOPPING CRITERIA
opt.maxfeval = 2000;        % maximum iteration count, in a sense
opt.tol = 1e-6;             % desired level of accuracy

%%% THE PROBLEM
opt.gamma = -4;             % the gamma thing
