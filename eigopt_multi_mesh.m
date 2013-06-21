%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hist, boxes] = eigopt_multi_mesh(funname, b0, b1, pars, set)
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
minmax = 0; % for now, let's just keep this simple
minmax = 1 - 2*boolean(minmax);

% the dimensions of the problem must be consistent
dim = length(b0);
dim2 = length(b1);
if (dim ~= dim2)
    error('lengths of b0 and b1 must be the same');
end

% generate the first 2^dim cells
xx0 = b0;
xx1 = b1;

% boxes contains the mesh cells
for j = 1:size(xx0,1),
    x0 = xx0(j,:);
    x1 = xx1(j,:);
    boxes(j) = initbox(funname,x0,x1,set.gamma,set.maxfeval,set.tol,pars,minmax);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If ever boxes(i).LB > min(boxes.UB), then box <i> can be discarded.

startflg = 1;
plotlast = 0;
displast = 0;
histindx = 0;
nfevals = 0;
hist = repmat(struct('nfevals',NaN,'nverts',NaN,'ncells', ...
    NaN,'LB',NaN,'UB',NaN,'z',NaN,'f',NaN, 'err', NaN), ...
    set.maxfeval, 1);

while startflg || ( abs(hist(histindx).UB - hist(histindx).LB) > set.tol ...
                    && hist(histindx).nfevals <= set.maxfeval ),
    startflg = 0;
    histindx = histindx + 1;

    % SPLIT if we are down to too few cells OR if a cell is past its prime
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if set.splitting,
        actives = find([boxes.active] == 1);
        splitlist = find([boxes(actives).iternum] > set.cellsmitosisage);
        splitlist = actives(splitlist);
        if length(actives) <= set.endangeredlimit,
            splitlist = actives;
        end
        
        for j = splitlist,
            boxes(j).active = -1;
            [newx0, newx1] = mitosis(boxes(j).lb, boxes(j).ub);
            for k = 1:2^dim,
                boxes(length(boxes)+1) = initbox(funname, ...
                    newx0(k,:),newx1(k,:),set.gamma,set.maxfeval,set.tol,pars,minmax);
            end
        end
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
    plotnow = nfevals - plotlast > set.plotfreq;
    dispnow = nfevals - displast > set.dispfreq;
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
        fprintf('nfevals:%d  LB:%.8f  UB:%.8f  #verts:%d  #cells:%d\n', ...
            hist(histindx).nfevals, hist(histindx).LB, hist(histindx).UB, ...
            hist(histindx).nverts, hist(histindx).ncells);
    end
end

% chop off the NaNs
indx = find(isnan([hist.UB]),1,'first');
hist(indx:end) = [];

% display final state
hist(end)

return
