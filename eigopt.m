%------------------------------------------------------------------%
% eigopt.m
% 
% An algorithm for minimizing an N-dimensional function.
%------------------------------------------------------------------%
function [fmin, xmin, hist, boxes] = eigopt(func, bounds, opts, varargin)

%-- run a demo function if no arguments are passed
if nargin == 0,
    func = 'fdist_uncont';
    bounds = [-2 6; -4 4];
    rng(1);
    pars.A = full(gallery('poisson',15));
    pars.B = randn(225,100);
    varargin = {pars};
end

%-- if the '-gethandles' flag is passed, return all the subfunction handles.
%-- this is for unit testing on the subfunctions.
if nargin == 1,
    if lower(func) == '-gethandles',
        fmin = {    @OPTbreadthsearch;  @OPTdepthsearch;  @OPTsetopts; ...
                    @OPTinitbox;  @OPTstep;  @OPTevalq;  @OPTisboundary; ...
                    @OPTmitosis;  @OPTheapsort;  @OPTheapinsert; ...
                    @OPTheapremove;  @OPTheapupdate;  @OPTplotgraph };
        return
    else
        error('Please specify bounds for the function.')
    end
end

%-- set user options
if ~exist('opts') || ~isstruct(opts), opts = []; end
defaults.plotfreq   = 40;       %-- frequency of plotting (if 2-D)
defaults.dispfreq   = 40;       %-- frequency of verbose output of status
defaults.splitting  = 1;        %-- recursively split into subboxes?
defaults.searchtype = 0;        %-- 0 == breadth-first; 1 == depth-first
defaults.maxquads   = 40;       %-- maximum number of quadratics per box
defaults.maxfeval   = 4000;     %-- maximum number of function evaluations
defaults.tol        = 1e-3;     %-- desired level of accuracy
defaults.gamma      = -1;       %-- the gamma thing
opts = OPTsetopts(opts, defaults);

%-- search for the solution via the specified method
if opts.searchtype == 0,
    [fmin,xmin,hist,boxes] = OPTbreadthsearch(func, bounds, opts, varargin);
else
    error('Depth-first search doesn''t work yet.')
    [fmin,xmin,hist,boxes] = OPTdepthsearch(func, bounds, opts, varargin);
end

%-- the final state
hist(end)

return


%------------------------------------------------------------------%
% Sets user options.
%
% opts      = a struct of options
% defaults  = a similar struct of default options
%
% Any options not defined in `opts` are set to the value in `defaults`.
%------------------------------------------------------------------%
function opts = OPTsetopts(opts, defaults)
ff = fieldnames(defaults);
for k = 1:length(ff),
    f = ff{k};
    if ~isfield(opts, f),
        opts.(f) = defaults.(f);
    end
end
return


%------------------------------------------------------------------%
% Breadth-first search for a global extremum of the function
%------------------------------------------------------------------%
function [fmin, xmin, hist, boxes] = OPTbreadthsearch(func, bounds, opts, varargin)

tic;                        %-- start the clock
dim = size(bounds,1);       %-- dimension of the problem

%-- `boxes` is a struct array of mesh boxes
%-- here we first preallocate for speed
boxes(1)      = OPTinitbox( func, bounds, opts, varargin{:} );
boxes(1000,1) = structfun(@(x) [], boxes(1), 'Uniform', 0);
boxcnt = 1;

% The Main Loop
% If ever  boxes(i).LB > boxes(obi).UB,  then box `i` can be discarded.
%------------------------------------------------------------------%
plotlast = 0;   %-- keep track of when to plot
displast = 0;   %-- keep track of when to verbose output
ii       = 0;   %-- keep track of the history

%-- save the history for analysis
hist = repmat(  struct( 'nfevals', NaN,  'nboxes', NaN, ...
                        'nverts',  NaN,  'x',      NaN, ...
                        'UB',      NaN,  'LB',     NaN, ...
                        'err',     NaN,  'toc',    NaN  ...
                    ), opts.maxfeval, 1);

%-- some global variables
nfevals = 1;            %-- total number of function evaluations
actives = 1;            %-- indices of active boxes
obi     = 1;            %-- index of the currently most optimal box

%-- The Main Loop
while abs(boxes(obi).UB - boxes(obi).LB) > opts.tol && nfevals < opts.maxfeval,
    ii = ii + 1;

    % SPLIT if we are down to too few cells OR if a cell is past its prime
    %------------------------------------------------------------------%
    if opts.splitting,
        if length(actives) == 1,
            splitlist = 1;
        elseif length(actives) > 1,
            splitlist = find([boxes(actives).iternum] >= opts.maxquads);
        else
            error('Something has gone seriously wrong. There are no active boxes.')
        end
        
        for j = splitlist,
            ji = actives(j);
            boxes(ji).status = -1;
            [b0, b1] = OPTmitosis(boxes(ji).bounds);

            actives = [actives boxcnt+1:boxcnt+2^dim];
            for k = 1:2^dim,
                boxcnt = boxcnt + 1;
                boxbounds = [b0(:,k) b1(:,k)];
                boxes(boxcnt) = OPTinitbox(func, ...
                                    boxbounds, opts, varargin{:});
                nfevals = nfevals + 1;
            end
        end
        actives(splitlist) = [];
    end

    % REMOVE any cells whose LBs exceed the global UB
    %------------------------------------------------------------------%
    obi = actives( find([boxes(actives).UB] == min([boxes(actives).UB]),1) );
    remlist = find( [boxes(actives).LB] > boxes(obi).UB );
    for j = remlist, % there is apparently no vectorized way of doing this
        boxes(actives(j)).status = 0;
    end
    actives(remlist) = [];

    % ITERATE on any active cells
    %------------------------------------------------------------------%
    for j = actives,
        boxes(j) = OPTstep(boxes(j));
    end
    nfevals = nfevals + length(actives);

    % REMOVE any cells whose LBs exceed the global UB
    %------------------------------------------------------------------%
    obi = actives( find([boxes(actives).UB] == min([boxes(actives).UB]),1) );
    remlist = find( [boxes(actives).LB] > boxes(obi).UB );
    for j = remlist, % there is apparently no vectorized way of doing this
        boxes(actives(j)).status = 0;
    end
    actives(remlist) = [];

    % DECIDE whether to plot or display or not
    %------------------------------------------------------------------%
    plotnow = nfevals - plotlast > opts.plotfreq;
    dispnow = nfevals - displast > opts.dispfreq;
    if plotnow, plotlast = nfevals; end
    if dispnow, displast = nfevals; end

    % PLOT all this beautiful stuff (if in 2 dimensions)
    %------------------------------------------------------------------%
    if dim == 2 && plotnow,
        for j = actives,
            if j == actives(1),
                figure(1);
                clf; hold on;
                plotbounds = boxes(1).bounds';
                axis(plotbounds(:)');
            end
            OPTplotgraph(boxes(j));
        end
    end

    % SAVE some of the information
    %------------------------------------------------------------------%
    hist(ii) = struct(  'nfevals',  nfevals, ...
                        'nboxes',   length(actives), ...
                        'nverts',   sum([boxes(actives).heaplength]), ...
                        'x',        boxes(obi).xx(:,end), ...
                        'UB',       boxes(obi).UB, ...
                        'LB',       boxes(obi).LB, ...
                        'err',      boxes(obi).UB - boxes(obi).LB, ...
                        'toc',      toc ...
                        );

    % OUTPUT information at regular intervals
    %------------------------------------------------------------------%
    if dispnow,
        fprintf('nfevals:%d  acc:%.1e  UB:%.8f  heaplen:%d  #boxes:%d (%d)\n', ...
            hist(ii).nfevals, hist(ii).err, ...
            hist(ii).UB, ...
            hist(ii).nverts, hist(ii).nboxes, ...
            max([boxes(actives).iternum]));
    end
end

hist(ii+1:end) = [];
boxes(boxcnt+1:end) = [];
fmin = boxes(obi).UB;
xmin = boxes(obi).xx(:,end);

return


%------------------------------------------------------------------%
% Depth-first search
%------------------------------------------------------------------%
function [fmin, xmin, hist, boxes] = OPTdepthsearch(func, bounds, opts, varargin)

dim = size(bounds,1);           %-- dimension of the problem

%-- `boxes` is a struct array of mesh boxes
boxes(1) = OPTinitbox(func, bounds, opts);

return


%------------------------------------------------------------------%
% initialize a new mesh box
%------------------------------------------------------------------%
function box = OPTinitbox(func, bounds, opts, varargin)

dim = size(bounds,1);
pars = varargin{:}{:};

box = struct( 'status', 1 ...
            ...%, 'f', @(x) feval(func, x, varargin{:}) ...
            , 'f', @(x) feval(func, x, pars) ...
            , 'bounds', bounds ...
            , 'gamma', opts.gamma ...
            , 'maxfeval', opts.maxfeval ...
            , 'tol', opts.tol ...
            , 'iternum', 1 ...
            , 'xx', [mean(bounds, 2) nan(dim,opts.maxfeval-1)] ...
            , 'fmin', nan(1, opts.maxfeval) ...
            , 'gmin', nan(dim, opts.maxfeval) ...
            , 'LB', -Inf ...
            , 'UB', Inf ...
            );

%-- evaluate the function and its gradient at x0
[box.fmin(1), box.gmin(:,1)] = box.f(box.xx(:,1));

% intialize vertices of the new box
%------------------------------------------------------------------%

%-- normalized boundary vertex coordinates
dim2 = 2^dim;
ee = ones(dim2,1);
coords = de2bi(0:dim2-1,dim);
coordscell = mat2cell( box.bounds(ee*(1:dim)+3*coords)', dim, ee );
[box.vertices(1:dim2).coor] = coordscell{:};

%-- boundary vertex index sets
indexcell = mat2cell( [ee*(-2*dim+1:2:-1)-fliplr(coords) ...
                        ee], ee );
[box.vertices(1:dim2).index] = indexcell{:};

%-- boundary vertex adjacencies
pows = 2.^(0:dim-1);
adj = bsxfun(@bitxor, [0:dim2-1]', pows)+1;     %-- generates hypercube adjacencies
sbadj = sparse(repmat(1:dim2, 1, dim), ...      %-- converts adj into a sparse
                adj(:), true, dim2, dim2);      %   logical matrix of adjacencies
sbadjcell = mat2cell(sbadj, ee);                %-- convert it to a cell array
[box.vertices(1:dim2).adj] = sbadjcell{:};      %-- finally, set all vertex adjacencies

%-- generate the quadratic model values at vertices
qq = OPTevalq([box.vertices.coor], box.xx(:,box.iternum), ...
               box.fmin(box.iternum), box.gmin(:,box.iternum), box.gamma);
qq = mat2cell(qq, 1, ee);
[box.vertices(1:dim2).qval] = qq{:};

%-- set up the heap for this box
box.heap = OPTheapsort(box.vertices);
box.heaplength = 2^dim;

%-- these are the initial bounds given by this box
box.LB = box.vertices(box.heap(1)).qval;
box.UB = box.fmin(1);

return


%------------------------------------------------------------------%
% iterate once on a mesh box
%------------------------------------------------------------------%
function box = OPTstep(box)

% INITIALIZATION
%------------------------------------------------------------------%

dim = size(box.bounds, 1);
box.iternum = box.iternum + 1;

%-- center for the new quadratic model
box.xx(:,box.iternum) = box.vertices(box.heap(1)).coor;

%-- function value and gradient at the center of the new model
[box.fmin(box.iternum),box.gmin(:,box.iternum)] = box.f(box.xx(:,box.iternum));

if (box.fmin(box.iternum) < box.UB)
    box.UB = box.fmin(box.iternum);
end

% (1) DETERMINE DEAD VERTICES
% some of these are interior dead
% some are dead on the boundary
%------------------------------------------------------------------%
stack(1) = box.heap(1);
stackl = 1;

notboundaryl = 0;
notboundarylist(1) = -1;
boundaryl = 0;
boundarylist(1) = -1;

while (stackl > 0),
    vertex = stack(stackl);
    stackl = stackl - 1;
    boundary = 0;

    %-- EXPAND TO THE ADJACENT VERTICES
    for adj = find(box.vertices(vertex).adj),

        %-- if the boolean is false, then `adj` is already inspected
        %   or to be inspected soon (in the stack)
        if ~any(adj == [boundarylist(1:boundaryl) ...
                notboundarylist(1:notboundaryl) stack(1:stackl)]),

            qnew = OPTevalq(box.vertices(adj).coor, box.xx(:,box.iternum), ...
                    box.fmin(box.iternum), box.gmin(:,box.iternum), box.gamma);

            if (qnew > box.vertices(adj).qval),
                %-- ADJACENT IS DEAD, ADD ADJACENT TO THE STACK
                stackl = stackl + 1;
                stack(stackl) = adj;
            else
                %-- ADJACENT IS ALIVE, SO VERTEX IS ON THE BOUNDARY OF THE SET OF DEAD VERTICES
                boundary = 1;
            end
        end
    end

    if  boundary == 0,
        %-- DEAD INTERIOR VERTICES
        notboundaryl = notboundaryl + 1;
        notboundarylist(notboundaryl) = vertex;
    else
        %-- DEAD BOUNDARY VERTICES
        boundaryl = boundaryl + 1;
        boundarylist(boundaryl) = vertex;
    end
end

deadboundaryl = 0;

% if (mod(box.iternum,10) == 0)
%     display(notboundaryl + boundaryl);
% end

% (2) REMOVE INTERIOR DEAD VERTICES
% Easy, no need to modify the alive
% vertices or introduce new vertices
%------------------------------------------------------------------%
for j = 1:notboundaryl

    if OPTisboundary(box.vertices,notboundarylist(j));
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist(deadboundaryl) = notboundarylist(j);
    else
        deadisboundaryvertex = 0;
    end

    if (deadisboundaryvertex == 0)
        % Dead vertex is not a boundary vertex, remove it
        %------------------------------------------------------------------%
        box.heap = OPTheapremove(box.heap,box.heaplength,box.vertices,notboundarylist(j));
        box.vertices(box.heap(box.heaplength)).qval = inf;

        box.heaplength = box.heaplength - 1;
    else
        % Dead vertex is a boundary vertex, cannot be removed
        % (i) but need to update the function value
        %------------------------------------------------------------------%
        deadboundary = notboundarylist(j);
        box.vertices(deadboundary).qval = OPTevalq(box.vertices(deadboundary).coor, ...
                                box.xx(:,box.iternum), box.fmin(box.iternum), ...
                                box.gmin(:,box.iternum),box.gamma);

        % (ii) and the index set
        %------------------------------------------------------------------%
        for k = 1:dim+1
            if ( box.vertices(deadboundary).index(k) > 0 )
                box.vertices(deadboundary).index(k) = box.iternum;
                break;
            end
        end

        box.vertices(deadboundary).adjnum = 0;

        % (iii) finally update the heap
        %------------------------------------------------------------------%
        box.heap = OPTheapupdate(box.heap,box.heaplength,box.vertices,deadboundary);
    end
end


% (3) REMOVE DEAD VERTICES ON THE BOUNDARY
% Subtle, (i) introduce a new vertex in between
% any pair of dead and alive vertices (ii) modify
% the adjacency of an alive vertex next to a dead
% vertex
%------------------------------------------------------------------%

verticesl = length(box.vertices);
oldlength = verticesl;
alldead = union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl));

%-- GO THROUGH THE LIST OF DEAD VERTICES ON THE BOUNDARY
for j = 1:boundaryl
    dead = boundarylist(j);
    qnewdead = OPTevalq(box.vertices(dead).coor, box.xx(:,box.iternum), ...
                     box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

    if OPTisboundary(box.vertices,dead);
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist(deadboundaryl) = dead;
    else
        deadisboundaryvertex = 0;
    end

    %-- Determine the adjacent alive vertices
    %-- FIXME:  try logical arrays instead:
    %           tic, nnz(setdiff(xx,intersect(xx,yy))), toc
    %           tic, nnz(sx & xor(sx,sy)), toc

    alivelist = setdiff(box.vertices(dead).adj, ...
                        intersect(box.vertices(dead).adj, alldead));
    alivelength = length(alivelist);

    for l = 1:alivelength
        alive = alivelist(l);

        %-- THIS ADJACENT VERTEX IS ALIVE, MUST CREATE A NEW VERTEX
        keyboard
        qnewalive = OPTevalq(box.vertices(alive).coor, box.xx(:,box.iternum), ...
                          box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

                % I - FORM A NEW VERTEX BETWEEN THE DEAD and ALIVE
                %------------------------------------------------------------------%
                verticesl = verticesl + 1;

                %(i) its coordinates
                %------------------------------------------------------------------%
                alpha = (qnewalive - box.vertices(alive).qval)/ ...
                    ((qnewalive - box.vertices(alive).qval) - (qnewdead - box.vertices(dead).qval));
                box.vertices(verticesl).coor = alpha*box.vertices(dead).coor ...
                                            + (1-alpha)*box.vertices(alive).coor;

                %(ii) its index
                %------------------------------------------------------------------%
                box.vertices(verticesl).index = intersect(box.vertices(alive).index, ...
                                                          box.vertices(dead).index);
                box.vertices(verticesl).index(dim+1) = box.iternum;

                %(iii) its function value
                %------------------------------------------------------------------%
                box.vertices(verticesl).qval = OPTevalq(box.vertices(verticesl).coor, ...
                                                     box.xx(:,box.iternum), ...
                                                     box.fmin(box.iternum), ...
                                                     box.gmin(:,box.iternum), box.gamma);
                for k = 1:box.iternum-1
                    qval = OPTevalq(box.vertices(verticesl).coor, ...
                                box.xx(:,k), box.fmin(k), box.gmin(:,k), box.gamma);
                    if (qval > box.vertices(verticesl).qval)
                        box.vertices(verticesl).qval = qval;
                    end
                end

                %(iv) initialize its adjacency list
                %------------------------------------------------------------------%
                box.vertices(verticesl).adjnum = 1;
                box.vertices(verticesl).adj(1) = alive;

                % II - ADD THE NEW VERTEX TO THE HEAP
                %------------------------------------------------------------------%
                box.heap = OPTheapinsert(box.heap,box.heaplength,box.vertices,verticesl);
                box.heaplength = box.heaplength + 1;

                % III - UPDATE THE ADJACENCY LIST OF THE ALIVE VERTEX
                %------------------------------------------------------------------%
                k = 1;
                while (box.vertices(alive).adj(k) ~= dead)
                    k = k+1;
                end
                box.vertices(alive).adj(k) = verticesl;
    end

    % IV - REMOVE THE DEAD VERTEX
    %------------------------------------------------------------------%
    if (deadisboundaryvertex == 0)
        %-- dead vertex is not a boundary vertex 
        box.heap = OPTheapremove(box.heap,box.heaplength,box.vertices,dead);
        box.vertices(dead).qval = inf;

        box.heaplength = box.heaplength - 1;
    else
        % dead vertex is on the boundary
        % (i) need to update the function value
        %------------------------------------------------------------------%
        box.vertices(dead).qval = qnewdead;

        % (ii) and the index set
        %------------------------------------------------------------------%
        for k = 1:dim+1
            if ( box.vertices(dead).index(k) > 0 )
                box.vertices(dead).index(k) = box.iternum;
                break;
            end
        end

        % (iii) and the adjacency
        %------------------------------------------------------------------%
        box.vertices(dead).adjnum = 0;

        % (iv) finally update the heap
        %------------------------------------------------------------------%
        box.heap = OPTheapupdate(box.heap,box.heaplength,box.vertices,dead);
    end
end

% (4) FORM THE CONNECTIONS IN THE NEW POLYTOPE
%------------------------------------------------------------------%
for j = oldlength+1:verticesl
    for k = j+1:verticesl
        if (length(intersect(box.vertices(j).index,box.vertices(k).index)) == dim)
            %-- newly added jth and kth vertices are adjacent
            adjnum = box.vertices(j).adjnum;
            box.vertices(j).adj(adjnum+1) = k;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;

            adjnum = box.vertices(k).adjnum;
            box.vertices(k).adj(adjnum+1) = j;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;    
        end
    end
end

for j = oldlength+1:verticesl
    for k = 1:deadboundaryl
        if (length(intersect(box.vertices(j).index,box.vertices(deadboundarylist(k)).index)) == dim)
            %-- newly added jth and the kth boundary vertices are adjacent
            adjnum = box.vertices(j).adjnum;
            box.vertices(j).adj(adjnum+1) = deadboundarylist(k);
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;
    
            adjnum = box.vertices(deadboundarylist(k)).adjnum;
            box.vertices(deadboundarylist(k)).adj(adjnum+1) = j;
            box.vertices(deadboundarylist(k)).adjnum = box.vertices(deadboundarylist(k)).adjnum + 1;    
        end
    end
end

for j = 1:deadboundaryl
    for k = j+1:deadboundaryl
        if (length(intersect(box.vertices(deadboundarylist(j)).index,...
                             box.vertices(deadboundarylist(k)).index)) == dim)
        %-- dead jth and the kth boundary vertices are adjacent
            adjnum = box.vertices(deadboundarylist(j)).adjnum;
            box.vertices(deadboundarylist(j)).adj(adjnum+1) = deadboundarylist(k);
            box.vertices(deadboundarylist(j)).adjnum = box.vertices(deadboundarylist(j)).adjnum + 1;

            adjnum = box.vertices(deadboundarylist(k)).adjnum;
            box.vertices(deadboundarylist(k)).adj(adjnum+1) = deadboundarylist(j);
            box.vertices(deadboundarylist(k)).adjnum = box.vertices(deadboundarylist(k)).adjnum + 1;    
        end
    end
end

box.LB = box.vertices(box.heap(1)).qval;
return


%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function quadval = OPTevalq(x,xk,fk,gk,gam)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified version July 30, 2012)

% quadval = fk + gk'*(x-xk) + (gam/2)*(norm(x-xk))^2;
ee = ones(1,size(x,2));
quadval = fk + gk'*(x-xk*ee) ...
        + (gam/2)*( sqrt(sum(abs(x-xk*ee).^2,1)) ).^2;
return;

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function bool = OPTisboundary(vertices,vertnum)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified version July 30, 2012)
bool = ~(length(find(vertices(vertnum).index > 0)) > 1);
return;

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
% divides a box into 2^dim smaller cells of equal size, e.g.
% > b = [0 1; 0 1];
% > bnew = OPTmitosis(b)
% `bnew` will have 2^dim rows, each corresponding to a new box
function [x0, x1] = OPTmitosis(b)
dim = size(b,1);
weights = de2bi((0:2^dim-1)')/2;
dx = ones(2^dim,1)*(b(:,2)-b(:,1))';
x0 = ones(2^dim,1)*b(:,1)' + weights.*dx;
x1 = x0 + dx/2;
x0 = x0'; x1 = x1';
return


%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function heap = OPTheapsort(vertices)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)
heap(1) = 1;
l = length(vertices);
for j = 2:l
    heap = OPTheapinsert(heap,j-1,vertices,j);
end
return;

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function heap = OPTheapinsert(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)
l = heaplength;
l = l+1;
heap(l) = index;
while (floor(l/2) > 0) & vertices(heap(l)).qval < vertices(heap(floor(l/2))).qval
    temp = heap(floor(l/2));
    heap(floor(l/2)) = heap(l);
    heap(l) = temp;
    l = floor(l/2);
end
return

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function heap = OPTheapremove(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)
index = find(heap == index);
vertices(heap(index)).qval = inf;

temp = heap(index);
heap(index) = heap(heaplength);
heap(heaplength) = temp;
heaplength = heaplength - 1;
while ( (2*index <= heaplength) ...
    & (vertices(heap(index)).qval > min(vertices(heap(2*index+1)).qval, ...
                                        vertices(heap(2*index)).qval)) )

    if (vertices(heap(2*index+1)).qval < vertices(heap(2*index)).qval)
        temp = heap(2*index + 1);
        heap(2*index + 1) = heap(index);
        heap(index) = temp;
        index = 2*index+1;
    else
        temp = heap(2*index);
        heap(2*index) = heap(index);
        heap(index) = temp;
        index = 2*index;
    end
end
return

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function [heap, heaplength] = OPTheapupdate(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)

% first find the vertex on the heap
index = find(heap == index);
while ((2*index <= heaplength) ...
        & (vertices(heap(index)).qval ...
            > min(vertices(heap(min(2*index+1,heaplength))).qval, ...
                  vertices(heap(2*index)).qval)))
    if (vertices(heap(min(2*index+1,heaplength))).qval ...
                            < vertices(heap(2*index)).qval)
        temp = heap(min(2*index + 1,heaplength));
        heap(min(2*index + 1,heaplength)) = heap(index);
        heap(index) = temp;
        index = min(2*index+1,heaplength);
    else
        temp = heap(2*index);
        heap(2*index) = heap(index);
        heap(index) = temp;
        index = 2*index;
    end
end
return


%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function OPTplotgraph(box)
if box.status ~= 1, return, end

% vertices and their adjacendies
for j = 1:box.heaplength
    cvertex = box.heap(j);
    
    for k = 1:box.vertices(cvertex).adjnum;
        avertex = box.vertices(cvertex).adj(k);   
        plot([box.vertices(cvertex).coor(1) box.vertices(avertex).coor(1)], ...
             [box.vertices(cvertex).coor(2) box.vertices(avertex).coor(2)],'k-');
    end
    
    plot(box.vertices(cvertex).coor(1),box.vertices(cvertex).coor(2),'k*');
end

% function evalation points
for j = 1:size(box.xx,2)
    plot(box.xx(1,j), box.xx(2,j), 'b*');
end

% the outline of the box
plot([box.bounds(1,1)*[1 1] box.bounds(1,2)*[1 1] box.bounds(1,1)], ...
     [box.bounds(2,1) box.bounds(2,s2)*[1 1] box.bounds(2,1)*[1 1]], 'c--')

axis off

return
