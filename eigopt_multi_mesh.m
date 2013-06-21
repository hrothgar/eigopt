%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hist, boxes] = eigopt_multi_mesh(funname, b0, b1, pars, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%    boxes(j).active
%          1 = iterating
%          0 = inactive; bounds proved it irrelevant
%         -1 = inactive; split into smaller cells
%
% BUGS
%   only supports minimization right now, not maximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert minmax to a sign constant
opt.minmax = 0; % for now, let's just keep this simple
opt.minmax = 1 - 2*boolean(opt.minmax);

% the dimensions of the problem must be consistent
dim = length(b0);
dim2 = length(b1);
if (dim ~= dim2)
    error('lengths of b0 and b1 must be the same');
end

% <boxes> is a struct array of mesh cells
boxes(1) = initbox(funname, b0, b1, pars, opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE MAIN LOOP
% If ever boxes(i).LB > min([boxes.UB]), then box <i> can be discarded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startflg = 1;   % just a flag for the first iteration
plotlast = 0;   % keep track of when to plot
displast = 0;   % keep track of when to verbose output
histindx = 0;   % keep track of the history

% save the history for analysis
hist = repmat(struct('nfevals', NaN, 'nverts', NaN, 'ncells', ...
    NaN, 'LB', NaN, 'UB', NaN, 'z', NaN, 'f', NaN, 'err', NaN), ...
    opt.maxfeval, 1);

% MAIN LOOP BEGIN
while startflg || ( abs(hist(histindx).UB - hist(histindx).LB) > opt.tol ...
                    && hist(histindx).nfevals <= opt.maxfeval ),
    startflg = 0;
    histindx = histindx + 1;

    % SPLIT if we are down to too few cells OR if a cell is past its prime
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.splitting,
        actives = find([boxes.active] == 1);
        splitlist = find([boxes(actives).iternum] > opt.maxquadspercell);
        splitlist = actives(splitlist);
        if length(actives) <= opt.mincellcount,
            splitlist = actives;
        end
        
        for j = splitlist,
            boxes(j).active = -1;
            [newx0, newx1] = mitosis(boxes(j).lb, boxes(j).ub);
            for k = 1:2^dim,
                boxes(length(boxes)+1) = initbox(funname, ...
                                    newx0(k,:), newx1(k,:), pars, opt);
            end
        end
    end

    % REMOVE any cells whose LBs exceed the global UB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indx = actives(find([boxes(actives).LB] > min([boxes.UB])));
    for i = indx, % there is apparently no vectorized way of doing this
        boxes(i).active = 0;
    end

    % ITERATE on any active cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    actives = find([boxes.active] == 1);
    for j = actives,
        boxes(j) = step(boxes(j));
    end

    % REMOVE any cells whose LBs exceed the global UB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indx = actives(find([boxes(actives).LB] > min([boxes.UB])));
    for i = indx, % there is apparently no vectorized way of doing this
        boxes(i).active = 0;
    end

    % DECIDE whether to plot or display or not
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfevals = sum([boxes.nfevals]);
    actives = find([boxes.active] == 1);
    plotnow = nfevals - plotlast > opt.plotfreq;
    dispnow = nfevals - displast > opt.dispfreq;
    if plotnow, plotlast = nfevals; end
    if dispnow, displast = nfevals; end

    % PLOT all this beautiful stuff (if in 2 dimensions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aboxes = boxes(actives);
    hist(histindx).nfevals = nfevals;
    hist(histindx).nverts = sum([aboxes.heaplength]);
    hist(histindx).ncells = length(actives);

    indx = find([aboxes.UB] == min([aboxes.UB]), 1);
    hist(histindx).LB = aboxes(indx).LB;
    hist(histindx).UB = aboxes(indx).UB;
    hist(histindx).z = aboxes(indx).quad(:,end);
    hist(histindx).err = hist(histindx).UB - hist(histindx).LB;
    hist(histindx).f = hist(histindx).UB; % might as well

    % OUTPUT information at regular intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dispnow,
        fprintf('nfevals:%d  acc:%.1e  UB:%.8f  #verts:%d  #cells:%d (%d)\n', ...
            hist(histindx).nfevals, hist(histindx).err, ...
            hist(histindx).UB, ...
            hist(histindx).nverts, hist(histindx).ncells, ...
            max([boxes(actives).iternum]));
    end
end

% chop off the NaNs
indx = find(isnan([hist.UB]),1,'first');
hist(indx:end) = [];

% display final state
hist(end)

return
