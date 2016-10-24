function [xd, yd, delta, deltaCI] = shifthd_pbci(x,y,q,nboot,alpha,plotit)
% [xd, yd, delta, deltaCI] = shifthd_pbci(x,y,q,nboot,alpha,plotit)
% Computes a shift function for two independent groups using the
% Harrell-Davis quantile estimator in conjunction with a percentile
% bootstrap approach.
%
% INPUTS:
% - x & y are vectors; no missing values are allowed
% - q = quantiles to estimate - default = deciles 0.1:0.1:.9
% - nboot = number of bootstrap samples - default = 2000
% - alpha = expected long-run type I error rate - default = 0.05
% - plotit = 1 to get a figure; 0 otherwise by default
%
% OUTPUTS:
% - xd & yd = vectors of quantiles
% - delta = vector of differences between quantiles (x-y)
% - deltaCI = matrix quantiles x low/high bounds of the confidence
% intervals of the quantile differences
%
% Unlike shifthd:
% - the confidence intervals are not corrected for multiple comparisons
% (the R version provides corrected critical p values - here i've decided
% not to provide p values at all - the goal is to understand how
% distributions differ, not to make binary decisions)
% - the confidence intervals are calculated using a percentile bootstrap of
%   the quantiles, instead of a percentile bootstrap of the standard error of the quantiles
% - the quantiles to compare can be specified and are not limited to the
%   deciles
% - Tied values are allowed
%
% Unlike Rand Wilcox's qcomhd R function, no p value is returned.
% Extensive experience suggests human beings cannot be trusted with p
% values.
%
% Adaptation of Rand Wilcox's qcomhd R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, SHIFTHD, SHIFTDHD_PBCI

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% GAR 2016-06-16 - first version

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
bootdelta = zeros(Nq,nboot);

lo = round(nboot.*alpha./2); % get CI boundary indices
hi = nboot - lo;
lo = lo+1;

for qi = 1:Nq
    xd(qi) = hd(x,q(qi));
    yd(qi) = hd(y,q(qi));
    delta(qi) = xd(qi) - yd(qi);
    for b = 1:nboot
        bootdelta(qi,b) = hd(x(randi(Nx,1,Nx)),q(qi)) - hd(y(randi(Ny,1,Ny)),q(qi));
    end
end

bootdelta = sort(bootdelta,2); % sort in ascending order
deltaCI(:,1) = bootdelta(:,lo);
deltaCI(:,2) = bootdelta(:,hi);

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
