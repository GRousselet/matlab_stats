function diffall_asym_res = diffall_asym(x,y,q,nboot,alpha,plotit)
% diffall_asym_res = diffall_asym(x,y,q,nboot,alpha,plotit)
% Computes a difference asymmetry function for two independent groups.
% See details in:
% Wilcox, R.R., Erceg-Hurn, D.M., Clark, F. & Carlson, M. (2014)
% Comparing two independent groups via the lower and upper quantiles.
% J Stat Comput Sim, 84, 1543-1551.
%
% The difference asymmetry function provides perspective on the degree a distribution
% is symmetric about zero, by quantifying the sum of q and 1-q quantiles of the
% distribution of all pairwise differences between two independent groups x and y.
%  If the distributions of x and y are identical, then the
%  distribution of all pairwise differences is symmetric about zero.
% If the distribution is symmetric the function should be approximately a horizontal line.
% If in addition the median of the difference scores is zero, the horizontal line will
% intersect the y-axis at zero.
% Confidence intervals and p values are returned for each quantile sum.
% The FWE is controlled via Hochberg's method, which is used to determine critical
% p values based on the argument alpha.
%
% INPUTS:
% - x and y are two vectors of independent observations
% - q = quantiles to estimate - default = 0.05:0.05:.4 - must be <.5
% - nboot = number of bootstrap samples - default = 1000
% - alpha = expected long-run type I error rate - default = 0.05
% - plotit = 1 to get a figure; 0 otherwise by default
%
% OUTPUTS:
% saved in structure diff_asym_res with fields:
% - q = estimated quantiles
% - hd_q = Harrell-Davis estimate of quantile q
% - hd_1minusq = Harrell-Davis estimate of quantile 1-q
% - qsum = sum of hd_q and hd_1minusq - the main quantity of interest
% - qsumci = matrix quantiles x low/high bounds of the confidence
%           intervals of the quantile sums
% - pval_crit = critical p value, corrected for multiple comparisons
% - pval = p value
%
% Adaptation of Rand Wilcox's qwmwhd and cbmhd R functions,
% from Rallfun-v31.txt
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, DIFF_ASYM, SHIFTHD_PBCI

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% GAR 2016-10-07 - first version
% GAR 2016-12-06 - sample quantiles independently

if nargin<3 || isempty(q)
    q = .05:.05:.40;
end
if nargin<4 || isempty(nboot)
    nboot = 1000;
end
if nargin<5 || isempty(alpha)
    alpha = 0.05;
end
if nargin<6 || isempty(plotit)
    plotit = 0;
end

x = x(:);
Nx = numel(x);
y = y(:);
Ny = numel(y);
Nq = length(q);
hd_q = zeros(Nq,1);
hd_1minusq = zeros(Nq,1);
qsum_ci = zeros(Nq,2);
pval = zeros(Nq,1);

% all pairwise differences
yy = repmat(y,[1 length(x)])';
xx = repmat(x,[1 length(y)]);
alldiff = xx-yy;

bootsamples = zeros(nboot,Nq);
for qi = 1:Nq
  % independent bootstrap samples for each quantile
  boottable_x = randi(Nx,Nx,nboot);
  boottable_y = randi(Ny,Ny,nboot);
    for B = 1:nboot
        yy = repmat(y(boottable_y(:,B)),[1 length(x)])';
        xx = repmat(x(boottable_x(:,B)),[1 length(y)]);
        bootalldiff = xx-yy;
        bootsamples(B,qi) = hd(bootalldiff(:),q(qi)) + hd(bootalldiff(:),1-q(qi));
    end
end

for qi = 1:Nq

hd_q(qi) = hd(alldiff(:),q(qi));
hd_1minusq(qi) = hd(alldiff(:),1-q(qi));

pv = mean(bootsamples(:,qi)<0) + .5*mean(bootsamples(:,qi)==0);
pval(qi) = 2*min(pv,1-pv);
low = round((alpha/2)*nboot)+1;
up = nboot-low;
sbvec = sort(bootsamples(:,qi));
qsum_ci(qi,1) = sbvec(low);
qsum_ci(qi,2) = sbvec(up);

end

qsum = hd_q + hd_1minusq;

temp = sort(pval,1,'descend');
zvec = alpha./(1:Nq);
pval_crit = zvec;

diffall_asym_res.q = q;
diffall_asym_res.hd_q = hd_q;
diffall_asym_res.hd_1minusq = hd_1minusq;
diffall_asym_res.qsum = qsum;
diffall_asym_res.qsum_ci = qsum_ci;
diffall_asym_res.pval_crit = pval_crit;
diffall_asym_res.pval = pval;

if plotit == 1

    figure('Color','w','NumberTitle','off');hold on

    ext = 0.1*(max(q)-min(q));
    plot([min(q)-ext max(q)+ext],[0 0],'LineWidth',1,'Color',[.5 .5 .5],'LineStyle','--') % zero line
    for qi = 1:Nq
        plot([q(qi) q(qi)],[qsum_ci(qi,1) qsum_ci(qi,2)],'k','LineWidth',2)
    end
    plot(q,qsum,'ko-','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10,'LineWidth',1)
    set(gca,'FontSize',14,'XLim',[min(q)-ext max(q)+ext])
    ext = 0.1*max(abs(qsum_ci(:)));
    set(gca,'YLim',[-max(abs(qsum_ci(:)))-ext max(abs(qsum_ci(:)))+ext])
    box on
    xlabel('Quantiles','FontSize',16,'FontWeight','bold')
    ylabel('Quantile sum = q + 1-q','FontSize',16,'FontWeight','bold')

end
