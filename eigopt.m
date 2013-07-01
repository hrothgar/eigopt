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
defaults.plotfreq   = 1;      %-- frequency of plotting (if 2-D)
defaults.dispfreq   = 10;       %-- frequency of verbose output of status
defaults.splitting  = 0;        %-- recursively split into subboxes?
defaults.searchtype = 0;        %-- 0 == breadth-first; 1 == depth-first
defaults.maxquads   = 40;       %-- maximum number of quadratics per box
defaults.maxfeval   = 4000;     %-- maximum number of function evaluations
defaults.tol        = 1e-6;     %-- desired level of accuracy
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
                        'x',        boxes(obi).xx(:,boxes(obi).iternum), ...
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
xmin = boxes(obi).xx(:,boxes(obi).iternum);

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

bignum = 5e3;               % FIXME: really
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
            , 'nvertices', NaN ...
            );

%-- evaluate the function and its gradient at x0
[box.fmin(1), box.gmin(:,1)] = box.f(box.xx(:,1));

% intialize vertices of the new box
%------------------------------------------------------------------%
%-- preallocate for vertices
box.vertices(1:bignum) = repmat( struct( ...
    'coor', NaN, 'index', NaN, ...
    'adj', sparse([],[],true,1,bignum,dim), ...
    'qval', NaN, 'adjnum', NaN ...
    ), 1, bignum);

%-- normalized boundary vertex coordinates
dim2 = 2^dim;
ee = ones(dim2,1);
coords = de2bi(0:dim2-1, dim);
coordscell = mat2cell( box.bounds(ee*(1:dim)+dim*coords)', dim, ee );
[box.vertices(1:dim2).coor] = coordscell{:};

%-- boundary vertex index sets
indexcell = mat2cell( [ee*(-2*dim+1:2:-1)-fliplr(coords) ee], ee );
[box.vertices(1:dim2).index] = indexcell{:};

%-- boundary vertex adjacencies
adj = bsxfun(@bitxor, [0:dim2-1]', 2.^(0:dim-1))+1; %-- generates hypercube adjacencies
sladj = sparse(repmat(1:dim2, 1, dim), ...          %-- converts adj into a sparse
                adj(:), true, dim2, bignum);        %   logical matrix of adjacencies
sladjcell = mat2cell(sladj, ee);                    %-- convert it to a cell array
[box.vertices(1:dim2).adj] = sladjcell{:};          %-- finally, set all vertex adjacencies

%-- generate the quadratic model values at vertices
qq = OPTevalq([box.vertices(1:dim2).coor], box.xx(:,box.iternum), ...
               box.fmin(box.iternum), box.gmin(:,box.iternum), box.gamma);
qq = mat2cell(qq, 1, ee);
[box.vertices(1:dim2).qval] = qq{:};

%-- adjacency numbers of each
adjnums = mat2cell(dim*ones(dim2,1),ones(dim2,1));
[box.vertices(1:dim2).adjnum] = adjnums{:};

box.nvertices = dim2;

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
box.UB = min(box.UB, box.fmin(box.iternum));

% (1) DETERMINE DEAD VERTICES
% some of these are interior dead
% some are dead on the boundary
%
% note: this could be faster maybe if written recursively
%------------------------------------------------------------------%
bignum = 5e3;               % FIXME: really
stack2 = sparse(1, box.heap(1), true, 1, bignum, 20);
stackl = 1;     % number of elements in stack

notboundaryl = 0;
notboundarylist2 = sparse([], [], true, 1, bignum, 20);

boundaryl = 0;
boundarylist2 = sparse([], [], true, 1, bignum, 20);

deadboundaryl = 0;
deadboundarylist2 = sparse([], [], true, 1, bignum, 20);

while (stackl > 0),
    vertex = find(stack2, 1, 'first');      %-- not optimal
    stack2(vertex) = false;
    stackl = stackl - 1;
    boundary = false;

    %-- EXPAND TO THE ADJACENT VERTICES
    for adj = find(box.vertices(vertex).adj),


        %-- if the boolean is false, then `adj` is already inspected
        %   or to be inspected soon (in the stack)
        if ~( boundarylist2(adj) | notboundarylist2(adj) | stack2(adj) ),

            qnew = OPTevalq(box.vertices(adj).coor, box.xx(:,box.iternum), ...
                    box.fmin(box.iternum), box.gmin(:,box.iternum), box.gamma);

            if (qnew > box.vertices(adj).qval),
                %-- ADJACENT IS DEAD, ADD ADJACENT TO THE STACK
                stackl = stackl + 1;
                stack2(adj) = true;
            else
                %-- ADJACENT IS ALIVE, SO VERTEX IS ON THE BOUNDARY OF THE SET OF DEAD VERTICES
                boundary = true;
            end
        end
    end

    if boundary,
        %-- DEAD BOUNDARY VERTICES
        boundaryl = boundaryl + 1;
        boundarylist2(vertex) = true;
    else
        %-- DEAD INTERIOR VERTICES
        notboundaryl = notboundaryl + 1;
        notboundarylist2(vertex) = true;
    end
end

% % number of dead vertices
% if (mod(box.iternum,10) == 0)
%     display(notboundaryl + boundaryl);
% end

% (2) REMOVE INTERIOR DEAD VERTICES
% Easy, no need to modify the alive
% vertices or introduce new vertices
%------------------------------------------------------------------%
for dead = find(notboundarylist2),
    if ~length(dead), break, end            %-- hack

    if OPTisboundary(box.vertices, dead);
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist2(dead) =  true;
    else
        deadisboundaryvertex = 0;
    end

    if deadisboundaryvertex,
        % Dead vertex is a boundary vertex, cannot be removed
        % (i) but need to update the function value
        %------------------------------------------------------------------%
        box.vertices(dead).qval = OPTevalq(box.vertices(dead).coor, ...
                                box.xx(:,box.iternum), box.fmin(box.iternum), ...
                                box.gmin(:,box.iternum),box.gamma);

        % (ii) and the index set
        %------------------------------------------------------------------%
        % FIXME yaaa
        for k = 1:dim+1,
            if ( box.vertices(dead).index(k) > 0 )
                box.vertices(dead).index(k) = box.iternum;
                break;
            end
        end
        box.vertices(dead).adjnum = 0;

        % (iii) finally update the heap
        %------------------------------------------------------------------%
    
        box.heap = OPTheapupdate(box.heap,box.heaplength,box.vertices,dead);
    else
        % Dead vertex is not a boundary vertex, remove it
        %------------------------------------------------------------------%
    
        box.heap = OPTheapremove(box.heap, box.heaplength, box.vertices, dead);

        %~ box.vertices(box.heap(box.heaplength)).qval = Inf;
        box.vertices(dead).qval = Inf;
        box.heaplength = box.heaplength - 1;
    end
end

% (3) REMOVE DEAD VERTICES ON THE BOUNDARY
% Subtle, (i) introduce a new vertex in between
% any pair of dead and alive vertices (ii) modify
% the adjacency of an alive vertex next to a dead
% vertex
%------------------------------------------------------------------%

oldlength = box.nvertices;
alldead2 = boundarylist2 | notboundarylist2;

%-- GO THROUGH THE LIST OF DEAD VERTICES ON THE BOUNDARY
for dead = find(boundarylist2),
    if ~length(dead), break, end           %-- hack

    qnewdead = OPTevalq(box.vertices(dead).coor, box.xx(:,box.iternum), ...
                     box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

    if OPTisboundary(box.vertices, dead);
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist2(dead) = true;
    else
        deadisboundaryvertex = 0;
    end

    %-- Determine the adjacent alive vertices
    alivelist2 = box.vertices(dead).adj & ~alldead2;

    for alive = find(alivelist2),
        if ~length(alive), 487, break, end           %-- hack

        %-- THIS ADJACENT VERTEX IS ALIVE, MUST CREATE A NEW VERTEX
        qnewalive = OPTevalq(box.vertices(alive).coor, box.xx(:,box.iternum), ...
                          box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

        % I - FORM A NEW VERTEX BETWEEN THE DEAD and ALIVE
        %------------------------------------------------------------------%
        box.nvertices = box.nvertices + 1;

        % (i) its coordinates
        %------------------------------------------------------------------%
        alpha = (qnewalive - box.vertices(alive).qval) / ...
            ( (qnewalive - box.vertices(alive).qval) - ...
                (qnewdead - box.vertices(dead).qval) );
        box.vertices(box.nvertices).coor = alpha*box.vertices(dead).coor ...
                                    + (1-alpha)*box.vertices(alive).coor;

        % (ii) its index
        %------------------------------------------------------------------%
        % FIXME, ya
        box.vertices(box.nvertices).index = intersect(box.vertices(alive).index, ...
                                                  box.vertices(dead).index);
        box.vertices(box.nvertices).index(dim+1) = box.iternum;

        % (iii) its function value
        %------------------------------------------------------------------%
        box.vertices(box.nvertices).qval = OPTevalq(box.vertices(box.nvertices).coor, ...
                                             box.xx(:,box.iternum), ...
                                             box.fmin(box.iternum), ...
                                             box.gmin(:,box.iternum), box.gamma);
        % if isnan(box.vertices(box.nvertices).qval), 505, keyboard, end

        % is this super inefficient? what is going on here?
        % FIXME WHAT IS THIS
        for k = 1:box.iternum-1,
            qval = OPTevalq(box.vertices(box.nvertices).coor, ...
                        box.xx(:,k), box.fmin(k), box.gmin(:,k), box.gamma);
            if (qval > box.vertices(box.nvertices).qval)
                box.vertices(box.nvertices).qval = qval;
            end
        end

        % (iv) initialize its adjacency list
        %------------------------------------------------------------------%
        box.vertices(box.nvertices).adjnum = 1;
        box.vertices(box.nvertices).adj(alive) = true;

        % II - ADD THE NEW VERTEX TO THE HEAP
        %------------------------------------------------------------------%
        box.heap = OPTheapinsert(box.heap,box.heaplength,box.vertices,box.nvertices);
        box.heaplength = box.heaplength + 1;

        % III - UPDATE THE ADJACENCY LIST OF THE ALIVE VERTEX
        %------------------------------------------------------------------%

        box.vertices(alive).adj(dead) = false;
        box.vertices(alive).adj(box.nvertices) = true;
    end

    % IV - REMOVE THE DEAD VERTEX
    %------------------------------------------------------------------%
    if deadisboundaryvertex,
        % dead vertex is on the boundary
        % (i) need to update the function value
        %------------------------------------------------------------------%
        box.vertices(dead).qval = qnewdead;

        % (ii) and the index set
        %------------------------------------------------------------------%
        % FIXME ya
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
    else
        %-- dead vertex is not a boundary vertex
        box.heap = OPTheapremove(box.heap,box.heaplength,box.vertices,dead);
        box.heaplength = box.heaplength - 1;
        box.vertices(dead).qval = Inf;
    end
end

% (4) FORM THE CONNECTIONS IN THE NEW POLYTOPE
%------------------------------------------------------------------%
for j = oldlength+1:box.nvertices
    for k = j+1:box.nvertices
        % FIXME eliminate inefficiencies
        if (length(intersect(box.vertices(j).index,box.vertices(k).index)) == dim)
            %-- newly added jth and kth vertices are adjacent
            box.vertices(j).adj(k) = true;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;

            box.vertices(k).adj(j) = true;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;    
        end
    end
end

for j = oldlength+1:box.nvertices
    for k = find(deadboundarylist2),
        if (length(intersect(box.vertices(j).index, box.vertices(k).index)) == dim)
            %-- newly added jth and the kth boundary vertices are adjacent
            box.vertices(j).adj(k) = true;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;
    
            box.vertices(k).adj(j) = true;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;
        end
    end
end

for j = find(deadboundarylist2),
    for k = find(deadboundarylist2(j:end)),
        if (length(intersect(box.vertices(j).index,...
                             box.vertices(k).index)) == dim)
        %-- dead jth and the kth boundary vertices are adjacent
            box.vertices(j).adj(k) = true;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;
            box.vertices(k).adj(j) = true;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;    
        end
    end
end

box.LB = box.vertices(box.heap(1)).qval;

% if any(isinf([box.vertices(box.heap).qval])), keyboard, end
% if any(isnan([box.vertices(box.heap).qval])), keyboard, end
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
l = find(isnan([vertices.adjnum]),1,'first') - 1;
for j = 2:l
    heap = OPTheapinsert(heap,j-1,vertices,j);
end
return;

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function heap = OPTheapinsert(heap,hl,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)
% hl = heaplength;
hl = hl+1;
heap(hl) = index;
while (floor(hl/2) > 0) & vertices(heap(hl)).qval < vertices(heap(floor(hl/2))).qval
    temp = heap(floor(hl/2));
    heap(floor(hl/2)) = heap(hl);
    heap(hl) = temp;
    hl = floor(hl/2);
end
return

%------------------------------------------------------------------%
% 
%------------------------------------------------------------------%
function heap = OPTheapremove(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)
ii = find(heap == index);
vertices(heap(ii)).qval = Inf;

temp = heap(ii);
heap(ii) = heap(heaplength);
heap(heaplength) = temp;
heaplength = heaplength - 1;
while ( (2*ii <= heaplength) ...
    & (vertices(heap(ii)).qval > min(vertices(heap(2*ii+1)).qval, ...
                                        vertices(heap(2*ii)).qval)) )

    if (vertices(heap(2*ii+1)).qval < vertices(heap(2*ii)).qval)
        temp = heap(2*ii + 1);
        heap(2*ii + 1) = heap(ii);
        heap(ii) = temp;
        ii = 2*ii+1;
    else
        temp = heap(2*ii);
        heap(2*ii) = heap(ii);
        heap(ii) = temp;
        ii = 2*ii;
    end
end
% heap = heap(1:end-1);
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
for j = 1:box.heaplength,
    cvertex = box.heap(j);

    for avertex = find(box.vertices(cvertex).adj);
        plot([box.vertices(cvertex).coor(1) box.vertices(avertex).coor(1)], ...
             [box.vertices(cvertex).coor(2) box.vertices(avertex).coor(2)],'k-');
    end

    plot(box.vertices(cvertex).coor(1),box.vertices(cvertex).coor(2),'k*');
end

% function evalation points
xx = [box.xx];
for j = 1:find(isnan(xx(1,:)), 1, 'first')-1,
    plot(xx(1,j), xx(2,j), 'b*');
end

% the outline of the box
plot([box.bounds(1,1)*[1 1] box.bounds(1,2)*[1 1] box.bounds(1,1)], ...
     [box.bounds(2,1) box.bounds(2,2)*[1 1] box.bounds(2,1)*[1 1]], 'c--')

axis off

return
