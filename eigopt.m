%------------------------------------------------------------------%
% eigopt.m
% 
% An algorithm for minimizing an N-dimensional function.
%------------------------------------------------------------------%
function [hist, boxes] = eigopt(func, bounds, opts, varargin)

if nargin < 3 || ~isstruct(opts), opts = []; end
defaults.plotfreq   = 20;       %-- frequency of plotting (if 2-D)
defaults.dispfreq   = 1;        %-- frequency of verbose output of status
defaults.splitting  = 1;        %-- recursively split into subboxes?
defaults.searchtype = 0;        %-- 0 == breadth-first; 1 == depth-first
defaults.maxquads   = 40;       %-- maximum number of quadratics per box
defaults.maxfeval   = 4000;     %-- maximum number of function evaluations
defaults.tol        = 1e-3;     %-- desired level of accuracy
defaults.gamma      = -4;       %-- the gamma thing
opts = setopts(opts, defaults);

%-- search for the solution via the specified method
if opts.searchtype == 0,
    [hist, boxes] = breadthsearch(func, bounds, opts, varargin);
else
    [hist, boxes] = depthsearch(func, bounds, opts, varargin);
end

return


%------------------------------------------------------------------%
% Breadth-first search
%------------------------------------------------------------------%
function [hist, boxes] = breadthsearch(func, bounds, opts, varargin)

%-- dimension of the problem
dim = size(bounds,1);

%-- `boxes` is a struct array of mesh boxes
%-- here we first preallocate for speed
boxes(1)       = initbox( func, bounds, opts, varargin{:} );
boxes(10000,2) = structfun(@(x) [], boxes(1), 'Uniform', 0);

% The Main Loop
% If ever  boxes(i).LB > gv.UB,  then box `i` can be discarded.
%------------------------------------------------------------------%
startflg = 1;   %-- just a flag to get things started
plotlast = 0;   %-- keep track of when to plot
displast = 0;   %-- keep track of when to verbose output
ii       = 0;   %-- keep track of the history

%-- save the history for analysis
hist = repmat(  struct( 'nfevals', NaN,  'nverts', NaN, ...
                        'nboxes',  NaN,  'LB',     NaN, ...
                        'UB',      NaN,  'z',      NaN, ...
                        'f',       NaN,  'err',    NaN  ...
                    ), opts.maxfeval, 1);

%-- some global variables
gv.UB = Inf;
gv.LB = -Inf;
gv.fevals = 1;

%-- The Main Loop
while startflg || ( abs(gv.UB - gv.LB) > opt.tol ...
                    && gv.nfevals < opt.maxfeval ),
    startflg = 0;
    ii = ii + 1;

    % SPLIT if we are down to too few cells OR if a cell is past its prime
    %------------------------------------------------------------------%
    if opt.splitting,
        actives = find([boxes.status] == 1);
        splitlist = actives(find([boxes(actives).iternum] >= opt.maxquads));
        if length(actives) <= 1, splitlist = actives; end
        
        for j = splitlist,
            boxes(j).status = -1;
            [newx0, newx1] = mitosis(boxes(j).lb, boxes(j).ub);
            for k = 1:2^dim,
                boxes(length(boxes)+1) = initbox(funname, ...
                                    [newx0(k,:) newx1(k,:)], opt, varargin{:});
            end
        end
    end

    % REMOVE any cells whose LBs exceed the global UB
    %------------------------------------------------------------------%
    actives = find([boxes.status] == 1);
    indx = actives(find([boxes(actives).LB] > min([boxes.UB])));
    for i = indx, % there is apparently no vectorized way of doing this
        boxes(i).status = 0;
    end

    % ITERATE on any active cells
    %------------------------------------------------------------------%
    actives = find([boxes.status] == 1);
    for j = actives,
        boxes(j) = step(boxes(j));
    end

    % REMOVE any cells whose LBs exceed the global UB
    %------------------------------------------------------------------%
    indx = actives(find([boxes(actives).LB] > min([boxes.UB])));
    for i = indx, % there is apparently no vectorized way of doing this
        boxes(i).status = 0;
    end

    % DECIDE whether to plot or display or not
    %------------------------------------------------------------------%
    nfevals = sum([boxes.nfevals]);
    actives = find([boxes.status] == 1);
    plotnow = nfevals - plotlast > opt.plotfreq;
    dispnow = nfevals - displast > opt.dispfreq;
    if plotnow, plotlast = nfevals; end
    if dispnow, displast = nfevals; end

    % PLOT all this beautiful stuff (if in 2 dimensions)
    %------------------------------------------------------------------%
    if dim == 2 && plotnow,
        for j = actives,
            if j == actives(1),
                figure(1);
                clf; hold on;
                axis([b0(1) b1(1) b0(2) b1(2)]);
            end
            plot_graph(boxes(j));
        end
    end

    % SAVE some of the information
    %------------------------------------------------------------------%
    aboxes = boxes(actives);
    hist(ii).nfevals = nfevals;
    hist(ii).nverts = sum([aboxes.heaplength]);
    hist(ii).ncells = length(actives);

    indx = find([aboxes.UB] == min([aboxes.UB]), 1);
    hist(ii).LB = aboxes(indx).LB;
    hist(ii).UB = aboxes(indx).UB;
    hist(ii).z = aboxes(indx).xx(:,end);
    hist(ii).err = hist(ii).UB - hist(ii).LB;
    hist(ii).f = hist(ii).UB; % might as well

    % OUTPUT information at regular intervals
    %------------------------------------------------------------------%
    if dispnow,
        fprintf('nfevals:%d  acc:%.1e  UB:%.8f  #verts:%d  #cells:%d (%d)\n', ...
            hist(ii).nfevals, hist(ii).err, ...
            hist(ii).UB, ...
            hist(ii).nverts, hist(ii).ncells, ...
            max([boxes(actives).iternum]));
    end
end

% chop off the NaNs
indx = find(isnan([hist.UB]),1,'first');
hist(indx:end) = [];

% display final state
hist(end)


return


%------------------------------------------------------------------%
% Depth-first search
%------------------------------------------------------------------%
function [hist, boxes] = depthsearch(func, bounds, opts, varargin)

dim = size(bounds,1);           %-- dimension of the problem

%-- `boxes` is a struct array of mesh boxes
boxes(1) = initbox(func, bounds, opts);

return


%------------------------------------------------------------------%
% Sets user options.
%
% opts      = a struct of options
% defaults  = a similar struct of default options
%
% Any options not defined in opts is set to the value in defaults.
%------------------------------------------------------------------%
function opts = setopts(opts, defaults)
ff = fieldnames(defaults);
for k = 1:length(ff),
    f = ff{k};
    if ~isfield(opts, f),
        opts.(f) = defaults.(f);
    end
end
return