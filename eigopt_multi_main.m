function [f, z, lbound, nfevals] = eigopt_multi_main(funname,lb,ub,gamma,itertol,tol,pars,fup,minmax)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim 
% (Modified version July 23, 2012)
%
% call: eigopt_multi_main(funname,lb,ub,gamma,itertol,tol,pars,fup,minmax)
% task:
%       this is the auxiliary routine called by eigopt_multi__mesh;
%       it applies the piece-wise quadratic model based algorithm
%       within the box [lb(1),ub(1)]x[lb(2),ub(2)]. Terminates
%       either 
%           (i) when the specified accuracy tolerance is satisfied, or
%           (ii) otherwise the number of quadratic functions in the
%           piece-wise model exceeds itertol
%
% Input:
%       funname (string) - name of the function evaluating 
%                  f(x) and gradient of f(x) at a given x
%       lb,ub (2x1 real) - the rectangle [lb(1),ub(1)]x[lb(2)xub(2)] must contain
%                          a global minimizer 
%       gamma (real)     - an upper bound for |phi''(alpha)| where
%                           phi(alpha) = f(x+alpha p) for all x and p
%       itertol (integer)   - maximum number of quadratic functions allowed
%       tol (real)   - the tolerance for the accuracy of the computed
%                  globally minimal value of f
%       pars (struct)    - parameters for funname
%       fup (real)      - a global upper bound for the globally minimal value of f
%       minmax (integer)    - minimize f(x) if minmax = 0
%                             otherwise maximize f(x)
% Output:
%       f (real)     - computed globally minimal value of f(x) 
%                  on [lb(1),ub(1)]x[lb(2),ub(2)]; differs from the exact 
%                  solution by no more than tol
%       z (2x1 real)  - computed global minimizer
%       nfevals (integer) - total number of function evaluations


% KEY DATA STRUCTURES
% vertices(j).
%   adjacency - array of indices to adjacent vertices
%   adjnum - no of adjacent integers
%   index - array of indices indicating active constraints
%   coor - vector of dim keeping the coordinates
%   quad - value of the largest quadratic function
%
% heap
%   keeps the indices of vertices sorted from the largest to the
%   smallest according to the value of the maximal quadratic function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a sign multiplier controlling whether we maximize or minimize the function
minmax = 1 - 2*boolean(minmax);

% the dimensions of the problem must be consistent
dim = length(lb);
dim2 = length(ub);
if (dim ~= dim2)
    error('lengths of lb and ub must be the same');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x0 = mid point of the box
quad = (lb(:)+ub(:))/2;

% evaluate the function and its gradient at x0
[fmin(1),gmin(:,1)] = feval(funname,quad(:,1),pars);
fmin(1) = minmax*fmin(1);
gmin(:,1) = minmax*gmin(:,1);

% construct the boundary vertices
for j = 1:2^dim
    k = j-1;
    l = 1;

    while (k ~= 0)
        if (mod(k,2) == 0)
            vertices(j).coor(l,1) = lb(l);
            vertices(j).index(l) = -(2*l - 1);
        else
            vertices(j).coor(l,1) = ub(l);
            vertices(j).index(l) = -2*l;
        end
        k = floor(k/2);
        l = l+1;
    end

    for k = l:dim
        vertices(j).coor(k,1) = lb(k);
        vertices(j).index(k) = -(2*k - 1);
    end

    vertices(j).adjnum = 0;
    vertices(j).quad = evalq(vertices(j).coor, quad(:,1), fmin(1), gmin(:,1),gamma);
    vertices(j).index(dim+1) = 1;
end

% boundary vertex adjacencies
for j = 1:2^dim
    for k = j+1:2^dim
        if (length(intersect(vertices(j).index,vertices(k).index)) == dim)
            adjnum = vertices(j).adjnum;
            vertices(j).adjacency(adjnum+1) = k;
            vertices(j).adjnum = vertices(j).adjnum + 1;

            adjnum = vertices(k).adjnum;
            vertices(k).adjacency(adjnum+1) = j;
            vertices(k).adjnum = vertices(k).adjnum + 1;
        end
    end
end

heap = heapsort(vertices);
heaplength = 2^dim;

lbound = vertices(heap(1)).quad;
ubound = fmin(1);

iternum = 1;

while ((ubound - lbound > tol) & (iternum <= itertol))

    % fprintf('iter:%d lowbound:%.8f upbound:%.8f heaplength:%d\n', ...
    %     iternum, minmax*lbound, minmax*ubound, heaplength);

    iternum = iternum + 1;

    % center for the new quadratic model
    quad(:,iternum) = vertices(heap(1)).coor;

    % function value and gradient at the center of the new model
    [fmin(iternum),gmin(:,iternum)] = feval(funname,quad(:,iternum),pars);

    fmin(iternum) = minmax*fmin(iternum);
    gmin(iternum) = minmax*gmin(iternum);

    if (fmin(iternum) < ubound)
        ubound = fmin(iternum);
        z = quad(:,iternum);
    end

    % (1) DETERMINE DEAD VERTICES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some of these are interior dead
    % some are dead on the boundary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stack(1) = heap(1);
    stackl = 1;

    notboundaryl = 0;
    notboundarylist(1) = -1;
    boundaryl = 0;
    boundarylist(1) = -1;

    while (stackl > 0)

        vertex = stack(stackl);
        stackl = stackl - 1;

        adjnum = vertices(vertex).adjnum;
        boundary  = 0;

        % EXPAND TO THE ADJACENT VERTICES
        for j = 1:adjnum

            adjacent = vertices(vertex).adjacency(j);

            if (ismember(adjacent, ...
                union(union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl)),stack(1:stackl)) ...
                    ) == 0)
                    % otherwise adjacent is already inspected or to be inspected soon (in the stack)

                qnew = evalq(vertices(adjacent).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

                if (qnew > vertices(adjacent).quad)
                    % ADJACENT IS DEAD, ADD ADJACENT TO THE STACK
                    stackl = stackl + 1;
                    stack(stackl) = adjacent;
                else
                    % ADJACENT IS ALIVE, SO VERTEX IS ON THE BOUNDARY OF THE SET OF DEAD VERTICES
                    boundary = 1;
                end
            end
        end

        if  boundary == 0
            % DEAD INTERIOR VERTICES
            notboundaryl = notboundaryl + 1;
            notboundarylist(notboundaryl) = vertex;
        else
            % DEAD BOUNDARY VERTICES
            boundaryl = boundaryl + 1;
            boundarylist(boundaryl) = vertex;
        end

    end

    deadboundaryl = 0;

    % (2) REMOVE INTERIOR DEAD VERTICES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Easy, no need to modify the alive
    % vertices or introduce new vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:notboundaryl

        if isboundary(vertices,notboundarylist(j));
            deadisboundaryvertex = 1;
            deadboundaryl = deadboundaryl + 1;
            deadboundarylist(deadboundaryl) = notboundarylist(j);
        else
            deadisboundaryvertex = 0;
        end

        if (deadisboundaryvertex == 0)
            % Dead vertex is not a boundary vertex, remove it
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            heap = heapremove(heap,heaplength,vertices,notboundarylist(j));
            vertices(heap(heaplength)).quad = inf;

            heaplength = heaplength - 1;
        else
            % Dead vertex is a boundary vertex, cannot be removed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (i) but need to update the function value
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            deadboundary = notboundarylist(j);
            vertices(deadboundary).quad = evalq(vertices(deadboundary).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

            %%%%%%%%%%%%%%%%%%%%%%%%
            % (ii) and the index set
            %%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:dim+1
                if ( vertices(deadboundary).index(k) > 0 )
                    vertices(deadboundary).index(k) = iternum;
                    break;
                end
            end

            vertices(deadboundary).adjnum = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (iii) finally update the heap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            heap = heapupdate(heap,heaplength,vertices,deadboundary);
        end
    end


    % (3) REMOVE DEAD VERTICES ON THE BOUNDARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subtle, (i) introduce a new vertex in between
    % any pair of dead and alive vertices (ii) modify
    % the adjacency of an alive vertex next to a dead
    % vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    verticesl = length(vertices);
    oldlength = verticesl;
    alldead = union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl));

    % GO THROUGH THE LIST OF DEAD VERTICES ON THE BOUNDARY
    for j = 1:boundaryl

        dead = boundarylist(j);
        adjnum = vertices(dead).adjnum;

        qnewdead = evalq(vertices(dead).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

        if isboundary(vertices,dead);
            deadisboundaryvertex = 1;
            deadboundaryl = deadboundaryl + 1;
            deadboundarylist(deadboundaryl) = dead;
        else
            deadisboundaryvertex = 0;
        end

        % Determine the adjacent alive vertices
        alivelist = setdiff(vertices(dead).adjacency, intersect(vertices(dead).adjacency, alldead));
        alivelength = length(alivelist);

        for l = 1:alivelength
            alive = alivelist(l);

            % THIS ADJACENT VERTEX IS ALIVE, MUST CREATE A NEW VERTEX

            qnewalive = evalq(vertices(alive).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % I - FORM A NEW VERTEX BETWEEN THE DEAD and ALIVE
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    verticesl = verticesl + 1;
                    %%%%%%%%%%%%%%%%%%%%
                    %(i) its coordinates
                    %%%%%%%%%%%%%%%%%%%%
                    alpha = (qnewalive - vertices(alive).quad)/((qnewalive - vertices(alive).quad) - (qnewdead - vertices(dead).quad));
                    vertices(verticesl).coor = alpha*vertices(dead).coor  +  (1-alpha)*vertices(alive).coor;

                    %%%%%%%%%%%%%%%%%%%%
                    %(ii) its index
                    %%%%%%%%%%%%%%%%%%%%
                    vertices(verticesl).index = intersect(vertices(alive).index,vertices(dead).index);
                    vertices(verticesl).index(dim+1) = iternum;

                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    %(iii) its function value
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    vertices(verticesl).quad = evalq(vertices(verticesl).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum), gamma);
                    for k = 1:iternum-1

                        qval = evalq(vertices(verticesl).coor, quad(:,k), fmin(k), gmin(:,k), gamma);
                        if (qval > vertices(verticesl).quad)
                            vertices(verticesl).quad = qval;
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %(iv) initialize its adjacency list
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    vertices(verticesl).adjnum = 1;
                    vertices(verticesl).adjacency(1) = alive;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % II - ADD THE NEW VERTEX TO THE HEAP
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    heap = heapinsert(heap,heaplength,vertices,verticesl);
                    heaplength = heaplength + 1;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % III - UPDATE THE ADJACENCY LIST OF THE ALIVE VERTEX
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    k = 1;
                    while (vertices(alive).adjacency(k) ~= dead)
                        k = k+1;
                    end
                    vertices(alive).adjacency(k) = verticesl;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IV - REMOVE THE DEAD VERTEX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (deadisboundaryvertex == 0)
            % dead vertex is not a boundary vertex 
            heap = heapremove(heap,heaplength,vertices,dead);
            vertices(dead).quad = inf;

            heaplength = heaplength - 1;
        else
            % dead vertex is on the boundary
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (i) need to update the function value
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vertices(dead).quad = qnewdead;

            %%%%%%%%%%%%%%%%%%%%%%%%
            % (ii) and the index set
            %%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:dim+1
                if ( vertices(dead).index(k) > 0 )
                    vertices(dead).index(k) = iternum;
                    break;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (iii) and the adjacency
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vertices(dead).adjnum = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (iv) finally update the heap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            heap = heapupdate(heap,heaplength,vertices,dead);
        end
    end

    % (4) FORM THE CONNECTIONS IN THE NEW POLYTOPE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = oldlength+1:verticesl
        for k = j+1:verticesl
        
            if (length(intersect(vertices(j).index,vertices(k).index)) == dim)
                % newly added jth and kth vertices are adjacent
                adjnum = vertices(j).adjnum;
                vertices(j).adjacency(adjnum+1) = k;
                vertices(j).adjnum = vertices(j).adjnum + 1;

                adjnum = vertices(k).adjnum;
                vertices(k).adjacency(adjnum+1) = j;
                vertices(k).adjnum = vertices(k).adjnum + 1;    
            end
        end
    end

    for j = oldlength+1:verticesl
        for k = 1:deadboundaryl
            if (length(intersect(vertices(j).index,vertices(deadboundarylist(k)).index)) == dim)
                % newly added jth and the kth boundary vertices are adjacent
                adjnum = vertices(j).adjnum;
                vertices(j).adjacency(adjnum+1) = deadboundarylist(k);
                vertices(j).adjnum = vertices(j).adjnum + 1;
        
                adjnum = vertices(deadboundarylist(k)).adjnum;
                vertices(deadboundarylist(k)).adjacency(adjnum+1) = j;
                vertices(deadboundarylist(k)).adjnum = vertices(deadboundarylist(k)).adjnum + 1;    
            end
        end
    end

    for j = 1:deadboundaryl
        for k = j+1:deadboundaryl
            if (length(intersect(vertices(deadboundarylist(j)).index,vertices(deadboundarylist(k)).index)) == dim)
            % dead jth and the kth boundary vertices are adjacent
                adjnum = vertices(deadboundarylist(j)).adjnum;
                vertices(deadboundarylist(j)).adjacency(adjnum+1) = deadboundarylist(k);
                vertices(deadboundarylist(j)).adjnum = vertices(deadboundarylist(j)).adjnum + 1;

                adjnum = vertices(deadboundarylist(k)).adjnum;
                vertices(deadboundarylist(k)).adjacency(adjnum+1) = deadboundarylist(j);
                vertices(deadboundarylist(k)).adjnum = vertices(deadboundarylist(k)).adjnum + 1;    
            end
        end
    end

    lbound = vertices(heap(1)).quad;

    % if (mod(iternum,10) == 0)
    %     plot_graph_main(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,1);
    % end

end

nfevals = iternum;
f = minmax*ubound;

return
