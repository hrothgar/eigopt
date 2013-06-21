function bool = isboundary(vertices,vertnum)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified version July 30, 2012)
%

bool = ~(length(find(vertices(vertnum).index > 0)) > 1);

return;