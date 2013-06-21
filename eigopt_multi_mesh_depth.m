%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hist, boxes] = eigopt_multi_mesh(funname, b0, b1, pars, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%    boxes(j).status
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
[newx0, newx1] = mitosis(b0, b1);
for k = 1:2^dim,
    boxes(k) = initbox(funname, ...
                    newx0(k,:), newx1(k,:), pars, opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE MAIN LOOP
% If ever boxes(i).LB > min([boxes.UB]), then box <i> can be discarded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPLETE = 0;   % completion flag
plotlast = 0;   % keep track of when to plot
displast = 0;   % keep track of when to verbose output
histindx = 1;   % keep track of the history

% save the history for analysis
hist = repmat(struct('nfevals', NaN, 'nverts', NaN, 'ncells', ...
    NaN, 'LB', NaN, 'UB', NaN, 'z', NaN, 'f', NaN, 'err', NaN), ...
    opt.maxfeval, 1);

% MAIN LOOP BEGIN
while 1,

    if COMPLETE, break; end

    actives = find([boxes.status] == 1);
    [~,bindx] = sort([boxes(actives).UB]);
    bindx = actives(bindx);

    for j = bindx,

        % ITERATE on this cell (if it's still relevant)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gUB = min([boxes.UB]);
        for k = 1:opt.maxquadspercell-1,
            if boxes(j).status ~= 1,
                break;
            end

            boxes(j) = step(boxes(j));

            % REMOVE this cell if it is a baddie
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if boxes(j).LB > gUB,
                boxes(j).status = 0;
            end
        end

        % DECIDE whether to plot or display or not
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nfevals = sum([boxes.nfevals]);
        actives = find([boxes.status] == 1);
        plotnow = nfevals - plotlast > opt.plotfreq;
        dispnow = nfevals - displast > opt.dispfreq;
        if plotnow, plotlast = nfevals; end
        if dispnow, displast = nfevals; end

        % PLOT all this beautiful stuff (if in 2 dimensions)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if dim == 2 && plotnow,
            for k = actives,
                if k == actives(1),
                    figure(1);
                    clf; hold on;
                    axis([b0(1) b1(1) b0(2) b1(2)]);
                end
                plot_graph(boxes(k));
            end
        end

        % SPLIT this box (if it's still relevant)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if boxes(j).status == 1,
            boxes(j).status = -1;
            [newx0, newx1] = mitosis(boxes(j).lb, boxes(j).ub);
            for k = 1:2^dim,
                boxes(length(boxes)+1) = initbox(funname, ...
                                    newx0(k,:), newx1(k,:), pars, opt);
            end
        end

        % REMOVE any cells whose LBs exceed the global UB
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        actives = find([boxes.status] == 1);
        indx = actives(find([boxes(actives).LB] > min([boxes.UB])));
        for k = indx, % there is apparently no vectorized way of doing this
            boxes(k).status = 0;
        end

        % SAVE some of the information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nfevals = sum([boxes.nfevals]);
        actives = find([boxes.status] == 1);
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
        
        % actives = find([boxes.status] == 1);
        % disp(['gUB:' num2str(min([boxes.UB])) ...
        %     '  gLB:' num2str(max([boxes(actives).LB])) ]);


        % CHECK the stopping criteria
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ( abs(hist(histindx).UB - hist(histindx).LB) <= opt.tol ...
                    || hist(histindx).nfevals > opt.maxfeval ),
            COMPLETE = 1;
            break;
        end

        % if hist(histindx).nfevals > 500, keyboard, end

        histindx = histindx + 1;

    end

end

% chop off the NaNs
indx = find(isnan([hist.UB]),1,'first');
hist(indx:end) = [];

% display final state
hist(end)

return
