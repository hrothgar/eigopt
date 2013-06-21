
function plot_results(r)

figure;

subplot(2,2,1)
plot([r.nfevals], [r.nverts])
xlabel('# function evals')
ylabel('# vertices'), title('# vertices')
xlim([min([r.nfevals]) max([r.nfevals])])

subplot(2,2,2)
plot([r.nfevals], [r.ncells])
ylim([min([r.ncells])-1 max([r.ncells])+1])
xlabel('# function evals')
ylabel('# mesh cells'), title('# mesh cells')
xlim([min([r.nfevals]) max([r.nfevals])])

subplot(2,2,3)
semilogy([r.nfevals], [r.UB]-[r.LB])
xlabel('# function evals')
ylabel('UB-LB'), title('UB-LB')
xlim([min([r.nfevals]) max([r.nfevals])])

subplot(2,2,4)
plot([r.nfevals], [r.f])
xlabel('# function evals')
ylabel('f_{min} estimate'), title('f_{min} estimate')
xlim([min([r.nfevals]) max([r.nfevals])])

end
