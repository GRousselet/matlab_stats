function [xhd,CI] = decilespbci(x,nboot,plotCI)
% [xhd CI] = decilespbci(x,nboot,plotCI)
% DECILESPBCI computes 95% percentile bootstrap confidence intervals
% for the deciles of a distribution using the Harrell-Davis estimator.
% Unlike decilesci, decilespbci does not rely on an estimation of the standard error of hd.
%
% INPUTS:
% x = vector
% nboot = number of bootstrap samples (default = 2000)
% plotCI = set to 1 to output figure of deciles + their confidence
% intervals
%
% OUTPUTS:
% xhd = 9 deciles
% CI = 9 x 2 matrix of the deciles' confidence intervals
%
% see:
% Wilcox, R.R. (2012)
% Introduction to robust estimation and hypothesis testing
% Academic Press
% p.126-132
%
% Adaptation of Rand Wilcox's qcipb & deciles R functions,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, HDCI, HDPBCI, DECILESCI

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% GAR 2016-06-02 - first version
% GAR 2016-06-13 - improved figure

if nargin<2;nboot=2000;plotCI=0;end

n = numel(x);
alpha = 0.05;
list = randi(n,nboot,n); % use same bootstrap samples for all CIs

% CI boundaries
lo = round(nboot*(alpha/2));
hi = nboot - lo;

xhd = zeros(9,1);
CI = zeros(9,2);

for d = 1:9
    
    q = d/10;
    
    % compute decile
    xhd(d) = hd(x,q);
    
    % percentile bootstrap of hd
    xboot = zeros(nboot,1);
    for B = 1:nboot
        xboot(B) = hd(x(list(B,:)),q);
    end
    xboot = sort(xboot);
    
    % confidence intervals
    CI(d,1) = xboot(lo+1);
    CI(d,2) = xboot(hi);
    
end

if plotCI==1
    
    figure('Color','w','NumberTitle','off');hold on
    plot([0 10],[0 0],'LineWidth',1,'Color',[.5 .5 .5]) % zero line
    for d = 1:9
       plot([d d],[CI(d,1) CI(d,2)],'k','LineWidth',2) 
    end
    % mark median
%     plot(5,xhd(5),'ko-','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10,'LineWidth',2)
    v = axis;plot([5 5],[v(3) v(4)],'k:')
    plot(1:9,xhd,'ko-','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10,'LineWidth',1)
%     plot(1:9,CI(:,1),'k+','MarkerFaceColor','k')
%     plot(1:9,CI(:,2),'k+','MarkerFaceColor','k')
    set(gca,'FontSize',14,'XLim',[0 10],'XTick',1:9)
    box on

end
