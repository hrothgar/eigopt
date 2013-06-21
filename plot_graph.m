function plot_graph(box)
if box.status ~= 1, return, end

% vertices and their adjacendies
for j = 1:box.heaplength
    cvertex = box.heap(j);
    
    for k = 1:box.vertices(cvertex).adjnum;
        avertex = box.vertices(cvertex).adjacency(k);   
        plot([box.vertices(cvertex).coor(1) box.vertices(avertex).coor(1)], ...
             [box.vertices(cvertex).coor(2) box.vertices(avertex).coor(2)],'k-');
    end
    
    plot(box.vertices(cvertex).coor(1),box.vertices(cvertex).coor(2),'k*');
end

% function evalation points
for j = 1:size(box.quad,2)
    plot(box.quad(1,j), box.quad(2,j), 'b*');
end

% the outline of the box
plot([box.lb(1)*[1 1] box.ub(1)*[1 1] box.lb(1)], ...
     [box.lb(2) box.ub(2)*[1 1] box.lb(2)*[1 1]], 'c--')

axis off

return
