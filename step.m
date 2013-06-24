%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate on a mesh cell ("box") one time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function box = step(box)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set.dispfreq = 10;                  % verbose output frequency
% set.plotfreq = 20;                  % replot every <plotfreq> iterations

dim = box.dim;

box.iternum = box.iternum + 1;

% center for the new quadratic model
box.xx(:,box.iternum) = box.vertices(box.heap(1)).coor;

% function value and gradient at the center of the new model
[box.fmin(box.iternum),box.gmin(:,box.iternum)] = box.f(box.xx(:,box.iternum));

if (box.fmin(box.iternum) < box.UB)
    box.UB = box.fmin(box.iternum);
    z = box.xx(:,box.iternum);
end

% (1) DETERMINE DEAD VERTICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some of these are interior dead
% some are dead on the boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stack(1) = box.heap(1);
stackl = 1;

notboundaryl = 0;
notboundarylist(1) = -1;
boundaryl = 0;
boundarylist(1) = -1;

while (stackl > 0)

    vertex = stack(stackl);
    stackl = stackl - 1;

    adjnum = box.vertices(vertex).adjnum;
    boundary  = 0;

    % EXPAND TO THE ADJACENT VERTICES
    for j = 1:adjnum

        adjacent = box.vertices(vertex).adjacency(j);

        if (ismember(adjacent, ...
            union( union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl)), ...
                    stack(1:stackl))) == 0)
                % otherwise adjacent is already inspected or to be inspected soon (in the stack)
            
            qnew = evalq(box.vertices(adjacent).coor, box.xx(:,box.iternum), ...
                    box.fmin(box.iternum), box.gmin(:,box.iternum), box.gamma);

            if (qnew > box.vertices(adjacent).xx)
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

% if (mod(box.iternum,10) == 0)
%     display(notboundaryl + boundaryl);
% end

% (2) REMOVE INTERIOR DEAD VERTICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Easy, no need to modify the alive
% vertices or introduce new vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:notboundaryl

    if isboundary(box.vertices,notboundarylist(j));
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist(deadboundaryl) = notboundarylist(j);
    else
        deadisboundaryvertex = 0;
    end

    if (deadisboundaryvertex == 0)
        % Dead vertex is not a boundary vertex, remove it
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box.heap = heapremove(box.heap,box.heaplength,box.vertices,notboundarylist(j));
        box.vertices(box.heap(box.heaplength)).xx = inf;

        box.heaplength = box.heaplength - 1;
    else
        % Dead vertex is a boundary vertex, cannot be removed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (i) but need to update the function value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deadboundary = notboundarylist(j);
        box.vertices(deadboundary).xx = evalq(box.vertices(deadboundary).coor, ...
                                box.xx(:,box.iternum), box.fmin(box.iternum), ...
                                box.gmin(:,box.iternum),box.gamma);

        %%%%%%%%%%%%%%%%%%%%%%%%
        % (ii) and the index set
        %%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:dim+1
            if ( box.vertices(deadboundary).index(k) > 0 )
                box.vertices(deadboundary).index(k) = box.iternum;
                break;
            end
        end

        box.vertices(deadboundary).adjnum = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iii) finally update the heap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box.heap = heapupdate(box.heap,box.heaplength,box.vertices,deadboundary);
    end
end


% (3) REMOVE DEAD VERTICES ON THE BOUNDARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtle, (i) introduce a new vertex in between
% any pair of dead and alive vertices (ii) modify
% the adjacency of an alive vertex next to a dead
% vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verticesl = length(box.vertices);
oldlength = verticesl;
alldead = union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl));

% GO THROUGH THE LIST OF DEAD VERTICES ON THE BOUNDARY
for j = 1:boundaryl

    dead = boundarylist(j);
    adjnum = box.vertices(dead).adjnum;

    qnewdead = evalq(box.vertices(dead).coor, box.xx(:,box.iternum), ...
                     box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

    if isboundary(box.vertices,dead);
        deadisboundaryvertex = 1;
        deadboundaryl = deadboundaryl + 1;
        deadboundarylist(deadboundaryl) = dead;
    else
        deadisboundaryvertex = 0;
    end

    % Determine the adjacent alive vertices
    alivelist = setdiff(box.vertices(dead).adjacency, intersect(box.vertices(dead).adjacency, alldead));
    alivelength = length(alivelist);

    for l = 1:alivelength
        alive = alivelist(l);

        % THIS ADJACENT VERTEX IS ALIVE, MUST CREATE A NEW VERTEX

        qnewalive = evalq(box.vertices(alive).coor, box.xx(:,box.iternum), ...
                          box.fmin(box.iternum), box.gmin(:,box.iternum),box.gamma);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % I - FORM A NEW VERTEX BETWEEN THE DEAD and ALIVE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                verticesl = verticesl + 1;
                %%%%%%%%%%%%%%%%%%%%
                %(i) its coordinates
                %%%%%%%%%%%%%%%%%%%%
                alpha = (qnewalive - box.vertices(alive).xx)/ ...
                    ((qnewalive - box.vertices(alive).xx) - (qnewdead - box.vertices(dead).xx));
                box.vertices(verticesl).coor = alpha*box.vertices(dead).coor ...
                                            + (1-alpha)*box.vertices(alive).coor;

                %%%%%%%%%%%%%%%%%%%%
                %(ii) its index
                %%%%%%%%%%%%%%%%%%%%
                box.vertices(verticesl).index = intersect(box.vertices(alive).index, ...
                                                          box.vertices(dead).index);
                box.vertices(verticesl).index(dim+1) = box.iternum;

                %%%%%%%%%%%%%%%%%%%%%%%%%
                %(iii) its function value
                %%%%%%%%%%%%%%%%%%%%%%%%%
                box.vertices(verticesl).xx = evalq(box.vertices(verticesl).coor, ...
                                                     box.xx(:,box.iternum), ...
                                                     box.fmin(box.iternum), ...
                                                     box.gmin(:,box.iternum), box.gamma);
                for k = 1:box.iternum-1

                    qval = evalq(box.vertices(verticesl).coor, ...
                                box.xx(:,k), box.fmin(k), box.gmin(:,k), box.gamma);
                    if (qval > box.vertices(verticesl).xx)
                        box.vertices(verticesl).xx = qval;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %(iv) initialize its adjacency list
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                box.vertices(verticesl).adjnum = 1;
                box.vertices(verticesl).adjacency(1) = alive;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % II - ADD THE NEW VERTEX TO THE HEAP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                box.heap = heapinsert(box.heap,box.heaplength,box.vertices,verticesl);
                box.heaplength = box.heaplength + 1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % III - UPDATE THE ADJACENCY LIST OF THE ALIVE VERTEX
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                k = 1;
                while (box.vertices(alive).adjacency(k) ~= dead)
                    k = k+1;
                end
                box.vertices(alive).adjacency(k) = verticesl;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IV - REMOVE THE DEAD VERTEX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (deadisboundaryvertex == 0)
        % dead vertex is not a boundary vertex 
        box.heap = heapremove(box.heap,box.heaplength,box.vertices,dead);
        box.vertices(dead).xx = inf;

        box.heaplength = box.heaplength - 1;
    else
        % dead vertex is on the boundary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (i) need to update the function value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box.vertices(dead).xx = qnewdead;

        %%%%%%%%%%%%%%%%%%%%%%%%
        % (ii) and the index set
        %%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:dim+1
            if ( box.vertices(dead).index(k) > 0 )
                box.vertices(dead).index(k) = box.iternum;
                break;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iii) and the adjacency
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box.vertices(dead).adjnum = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iv) finally update the heap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box.heap = heapupdate(box.heap,box.heaplength,box.vertices,dead);
    end
end

% (4) FORM THE CONNECTIONS IN THE NEW POLYTOPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = oldlength+1:verticesl
    for k = j+1:verticesl
        if (length(intersect(box.vertices(j).index,box.vertices(k).index)) == dim)
            % newly added jth and kth vertices are adjacent
            adjnum = box.vertices(j).adjnum;
            box.vertices(j).adjacency(adjnum+1) = k;
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;

            adjnum = box.vertices(k).adjnum;
            box.vertices(k).adjacency(adjnum+1) = j;
            box.vertices(k).adjnum = box.vertices(k).adjnum + 1;    
        end
    end
end

for j = oldlength+1:verticesl
    for k = 1:deadboundaryl
        if (length(intersect(box.vertices(j).index,box.vertices(deadboundarylist(k)).index)) == dim)
            % newly added jth and the kth boundary vertices are adjacent
            adjnum = box.vertices(j).adjnum;
            box.vertices(j).adjacency(adjnum+1) = deadboundarylist(k);
            box.vertices(j).adjnum = box.vertices(j).adjnum + 1;
    
            adjnum = box.vertices(deadboundarylist(k)).adjnum;
            box.vertices(deadboundarylist(k)).adjacency(adjnum+1) = j;
            box.vertices(deadboundarylist(k)).adjnum = box.vertices(deadboundarylist(k)).adjnum + 1;    
        end
    end
end

for j = 1:deadboundaryl
    for k = j+1:deadboundaryl
        if (length(intersect(box.vertices(deadboundarylist(j)).index,...
                             box.vertices(deadboundarylist(k)).index)) == dim)
        % dead jth and the kth boundary vertices are adjacent
            adjnum = box.vertices(deadboundarylist(j)).adjnum;
            box.vertices(deadboundarylist(j)).adjacency(adjnum+1) = deadboundarylist(k);
            box.vertices(deadboundarylist(j)).adjnum = box.vertices(deadboundarylist(j)).adjnum + 1;

            adjnum = box.vertices(deadboundarylist(k)).adjnum;
            box.vertices(deadboundarylist(k)).adjacency(adjnum+1) = deadboundarylist(j);
            box.vertices(deadboundarylist(k)).adjnum = box.vertices(deadboundarylist(k)).adjnum + 1;    
        end
    end
end

box.LB = box.vertices(box.heap(1)).xx;

box.nfevals = box.iternum;

return
