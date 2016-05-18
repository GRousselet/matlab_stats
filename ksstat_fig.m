function ks = ksstat_fig(x,y)
% ks = ksstat_fig(x,y)
% Computes & illustrate the two-sample Kolmogorov-Smirnov statistic.
% x & y are two independent groups.
% Returns ks, the maximum of the absolute differences between the two empirical
% cumulative distribution functions (CDF).
%
% No p value is returned.
%
% Based on Rand Wilcox's ks & ecdf R functions
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also KSSTAT AKERD

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% remove NaNs & reformat
x=x(~isnan(x)); x=x(:); Nx = numel(x);
y=y(~isnan(y)); y=y(:); Ny = numel(y);

all = sort([x;y]); % pool and sort the observations
ecdfx = zeros(size(all));
ecdfy = zeros(size(all));
for e = 1:numel(all)
    ecdfx(e) = length(x(x<=all(e)))/Nx;
    ecdfy(e) = length(y(y<=all(e)))/Ny;
end
diff = abs( ecdfx - ecdfy );
ks = max( diff );
xmax = find(diff == ks);
xmax = xmax(1);

xkde = akerd(x);
ykde = akerd(y);

figure('Color','w','NumberTitle','off')

subplot(2,1,1);hold on % kernel density estimates
plot(sort(x),xkde,'Color',[1 0.5 0.2],'LineWidth',2)
plot(sort(y),ykde,'Color',[0 .5 0],'LineWidth',2)
xlabel('x','FontSize',16)
ylabel('density','FontSize',16)

subplot(2,1,2);hold on % KS statistics
title(sprintf('KS = %.2f',ks),'FontSize',18)
plot([all(xmax) all(xmax)],[0 1],'k:')
plot([all(xmax) all(xmax)],[ecdfx(xmax) ecdfy(xmax)],'k','LineWidth',3)
plot(all,ecdfx,'Color',[1 0.5 0.2],'LineWidth',2)
plot(all,ecdfy,'Color',[0 .5 0],'LineWidth',2)
xlabel('x','FontSize',16)
ylabel('cumulative distribution','FontSize',16)

for sub = 1:2
    subplot(2,1,sub)
    box on
    axis tight
    set(gca,'FontSize',14,'LineWidth',1,'Layer','Top')
    if sub == 1
       v=axis;
       set(gca,'YLim',[0 v(4)*1.05])
    end
end