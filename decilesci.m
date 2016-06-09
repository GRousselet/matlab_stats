function [xhd,CI] = decilesci(x,nboot,plotCI)
% [xhd CI] = decilesci(x,nboot,plotCI)
% DECILESCI computes 95% confidence intervals
% for the deciles of a distribution using
% the Harrell-Davis estimator and a percentile bootstrap
% estimate of its standard error.
%
% INPUTS:
% x = vector
% nboot = number of bootstrap samples (default = 100)
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
% Adaptation of Rand Wilcox's hdci R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, HDCI, HDPBCI, DECILESPBCI

% Copyright (C) 2009, 2016 Guillaume Rousselet - University of Glasgow
% GAR, University of Glasgow, Sep 2009 - first version
% GAR 2016-06-02 - updated help

if nargin<2;nboot=100;plotCI=0;end

n = numel(x);
if n<=10
    error('hdci: n<11, confidence intervals cannot be computed for less than 11 observations')
end

list = randi(n,nboot,n); % use same bootstrap samples for all CIs

xhd = zeros(9,1);
CI = zeros(9,2);

for d = 1:9
    
    q = d/10;
    
    % compute decile
    xhd(d) = hd(x,q); 
    
    % percentile bootstrap estimate of the standard error of hd
    xboot = zeros(nboot,1);
    for B = 1:nboot
        xboot(B) = hd(x(list(B,:)),q);
    end
    xd_bse = std(xboot,0); % normalize by (n-1)
    % We estimate the standard error of the test statistic by the
    % standard deviation of the bootstrap replications
    % Efron & Tibshinari 1993, chapter 6
    % Wilcox 2005, p.44-45
    
    % The constant c was determined so that the probability coverage of each confidence interval
    % is approximately 95% when sampling from normal and non-normal distributions.
    c = 1.96 + .5064 ./ (n.^.25);
    
    if q<=.2 || q>=.8
        if n <= 20
            c = -6.23./n+5.01;
        end
    end
    
    if q<=.1 || q>=.9
        if  n<=40
            c = 36.2./n+1.31;
        end
    end
      
    % confidence intervals
    CI(d,1) = xhd(d)-c.*xd_bse;
    CI(d,2) = xhd(d)+c.*xd_bse;
    
end

if plotCI==1
    figure;set(gcf,'Color','w');hold on
    plot(1:9,xhd,'ko',1:9,CI(:,1),'k+',1:9,CI(:,2),'k+')
    set(gca,'FontSize',14,'XLim',[0 10],'XTick',1:9)
    box on
end
