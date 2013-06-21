function quadval = plot_graph(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,fignum)


figure(fignum);
hold off;

for j = 1:heaplength
	
	cvertex = heap(j);
	
	for k = 1:vertices(cvertex).adjnum;
		avertex = vertices(cvertex).adjacency(k);	
		plot([vertices(cvertex).coor(1) vertices(avertex).coor(1)],[vertices(cvertex).coor(2) vertices(avertex).coor(2)],'k-');  

		if ((j == 1) & (k == 1))
			hold on;
		end
	end

	
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'k*');
end


for j = 1:boundaryl
	cvertex = boundarylist(j);
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'b*');
end

for j = 1:notboundaryl
	cvertex = notboundarylist(j);
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'bs');
end



return;