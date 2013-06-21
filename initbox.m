%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize a new mesh cell ("box")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function box = initbox(funname,lb,ub,pars,opt)

box = struct('status', 1, 'funname', funname, ...
    'f', @(x) feval(funname, x, pars), 'lb', lb, 'ub', ub, ...
    'gamma', opt.gamma, 'maxfeval', opt.maxfeval, 'tol', opt.tol, ...
    'pars', pars, 'minmax', opt.minmax, 'dim', length(lb), ...
    'nfevals', 1, 'fval', opt.minmax*Inf);

dim = box.dim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE GRAPH
% this is everything before the main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x0 = mid point of the box
box.quad = (lb(:)+ub(:))/2;

% evaluate the function and its gradient at x0

[box.fmin(1),box.gmin(:,1)] = box.f(box.quad(:,1));
box.fmin(1) = box.minmax*box.fmin(1);
box.gmin(:,1) = box.minmax*box.gmin(:,1);

% construct the boundary vertices
for j = 1:2^dim
    k = j-1;
    l = 1;

    while (k ~= 0)
        if (mod(k,2) == 0)
            box.vertices(j).coor(l,1) = box.lb(l);
            box.vertices(j).index(l) = -(2*l - 1);
        else
            box.vertices(j).coor(l,1) = ub(l);
            box.vertices(j).index(l) = -2*l;
        end
        k = floor(k/2);
        l = l+1;
    end

    for k = l:dim
        box.vertices(j).coor(k,1) = box.lb(k);
        box.vertices(j).index(k) = -(2*k - 1);
    end

    box.vertices(j).adjnum = 0;
    box.vertices(j).quad = evalq(box.vertices(j).coor, ...
            box.quad(:,1), box.fmin(1), box.gmin(:,1), box.gamma);
    box.vertices(j).index(dim+1) = 1;
end

% boundary vertex adjacencies
for j = 1:2^dim
    for k = j+1:2^dim
        if (length(intersect(box.vertices(j).index,box.vertices(k).index)) == dim)
            adjnum = box.vertices(j).adjnum;
            box.vertices(j).adjacency(adjnum+1) = k;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;

            adjnum = box.vertices(k).adjnum;
            box.vertices(k).adjacency(adjnum+1) = j;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;
        end
    end
end

box.heap = heapsort(box.vertices);
box.heaplength = 2^dim;

box.LB = box.vertices(box.heap(1)).quad;
box.UB = box.fmin(1);

box.iternum = 1;

return
