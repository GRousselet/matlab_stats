function [xhd,CI] = hdci(x,q,nboot)
% [xhd CI] = hdci(x,q,nboot)
% HDCI computes the 95% confidence interval
% for the qth quantile of a distribution using
% the Harrell-Davis estimator and a percentile bootstrap
% estimate of its standard error.
%
% INPUTS:
% x = vector
% q = quantile (default 0.5)
% nboot = number of bootstrap samples (default = 100)
%
% OUTPUTS:
% xhd = quantile
% CI = quantile's confidence interval
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
% See also HD, HDPBCI, DECILESPBCI, DECILESCI

% Copyright (C) 2009, 2016 Guillaume Rousselet - University of Glasgow
% GAR, University of Glasgow, Sep 2009 - first version
% GAR 2016-06-02 - updated help

if nargin<2;q=0.5;nboot=100;end

n = numel(x);
if n<=10
    error('hdci: n<11, confidence intervals cannot be computed for less than 11 observations')
end

% compute decile
xhd = hd(x,q);

% percentile bootstrap estimate of the standard error of hd
xboot = zeros(nboot,1);
list = randi(n,nboot,n);
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
CI(1) = xhd - c.*xd_bse;
CI(2) = xhd + c.*xd_bse;
