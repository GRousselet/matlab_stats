function [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(x,y,q,nboot,alpha,plotit)
% [xd, yd, delta, deltaCI] = shiftdhd_pbci(x,y,q,nboot,alpha,plotit)
% Computes a shift function for two dependent groups by comparing 
% the quantiles of the marginal distributions using the
% Harrell-Davis quantile estimator in conjunction with a percentile
% bootstrap approach.
%
% INPUTS:
% - x & y are vectors of the same length without missing values
% - q = quantiles to estimate - default = 0.1:0.1:.9 (deciles)
% - nboot = number of bootstrap samples - default = 2000
% - alpha = expected long-run type I error rate - default = 0.05
% - plotit = 1 to get a figure; 0 otherwise by default
%
% OUTPUTS:
% - xd & yd = vectors of quantiles
% - delta = vector of differences between quantiles (x-y)
% - deltaCI = matrix quantiles x low/high bounds of the confidence
% intervals of the quantile differences
% - pval = vector of p values
% - cpval = vector of critical p values based on Hochberg's method
% - sig = logical vector indicating if pval <= critical pval
%
% Unlike shiftdhd:
% - the confidence intervals are not corrected for multiple comparisons
% - the confidence intervals are calculated using a percentile bootstrap of
%   the quantiles, instead of a percentile bootstrap of the standard error of the quantiles
% - the quantiles to compare can be specified and are not limited to the
%   deciles
% - Tied values are allowed
% - alpha is not restricted to 0.05
%
% Adaptation of Rand Wilcox's Dqcomhd R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% Reference:
% Wilcox, R.R. & Erceg-Hurn, D.M. (2012)
% Comparing two dependent groups via quantiles.
% J Appl Stat, 39, 2655-2664.
%
% See also HD, SHIFTDHD, SHIFTHD_PBCI

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% GAR 2016-06-16 - first version
% GAR 2017-03-30 - add p values and CI correction

if nargin<3 || isempty(q)
    q = .1:.1:.9;
end
if nargin<4 || isempty(nboot)
    nboot = 2000;
end
if nargin<5 || isempty(alpha)
    alpha = 0.05;
end
if nargin<6 || isempty(plotit)
    plotit = 0;
end

Nq = numel(q);
xd = zeros(Nq,1);
yd = zeros(Nq,1);
delta = zeros(Nq,1);
deltaCI = zeros(Nq,2);
Nx = numel(x);
Ny = numel(y);
if Nx ~= Ny
   error('shiftdhd_pbci: x and y must be the same length') 
end
bootdelta = zeros(Nq,nboot);

lo = round(nboot.*alpha./2); % get CI boundary indices
hi = nboot - lo;
lo = lo+1;

for qi = 1:Nq
    xd(qi) = hd(x,q(qi));
    yd(qi) = hd(y,q(qi));
    delta(qi) = xd(qi) - yd(qi);
    % independent bootstrap samples for each quantile
    % this is necessary to control false positives
    for b = 1:nboot
        bootsample = randi(Nx,1,Nx); % sample pairs of observations
        bootdelta(qi,b) = hd(x(bootsample),q(qi)) - hd(y(bootsample),q(qi));
    end
end

% confidence intervals
bootdelta = sort(bootdelta,2); % sort in ascending order
deltaCI(:,1) = bootdelta(:,lo);
deltaCI(:,2) = bootdelta(:,hi);

% p values
pval = sum(bootdelta<0,2)./nboot + sum(bootdelta==0,2)./(2*nboot);
pval = 2*(min(pval,1-pval));

% critical p values
[~,I] = sort(pval,1,'descend');
zvec = alpha./(1:Nq);
cpval = zvec(I);

sig = pval' <= cpval;

if plotit == 1
    
    figure('Color','w','NumberTitle','off');hold on
    
    ext = 0.1*(max(xd)-min(xd));
    plot([min(xd)-ext max(xd)+ext],[0 0],'LineWidth',1,'Color',[.5 .5 .5]) % zero line
    for qi = 1:Nq
       plot([xd(qi) xd(qi)],[deltaCI(qi,1) deltaCI(qi,2)],'k','LineWidth',2) 
    end
    % mark median
    % v = axis;plot([xd(5) xd(5)],[v(3) v(4)],'k:')
    plot(xd,delta,'ko-','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10,'LineWidth',1)
    set(gca,'FontSize',14,'XLim',[min(xd)-ext max(xd)+ext])
    box on
    xlabel('Group 1 quantiles','FontSize',16)
    ylabel('Group 1 - group 2 quantiles','FontSize',16)
    
end