function [xhd,CI] = hdpbci(x,q,nboot,alpha)
% [xhd CI] = hdpbci(x,q,nboot)
% HDPBCI computes the 95% percentile bootstrap confidence interval
% for the qth quantile of a distribution using the Harrell-Davis estimator.
% Unlike hdci, hdpbci does not rely on an estimation of the standard error of hd.
%
% INPUTS:
% x = vector
% q = quantile (default 0.5)
% nboot = number of bootstrap samples (default = 2000)
% alpha = alpha level (default = 0.05)
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
% Adaptation of Rand Wilcox's qcipb R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, HDCI, DECILESPBCI, DECILESCI

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% GAR 2016-06-02 - first version

if ~exist('q', 'var') || isempty(q)
    q = 0.5;
end
if ~exist('nboot', 'var') || isempty(nboot)
    nboot = 2000;
end
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 0.05;
end

n = numel(x);

% compute quantile
xhd = hd(x,q);

% percentile bootstrap of hd
xboot = zeros(nboot,1);
list = randi(n,nboot,n);
for B = 1:nboot
    xboot(B) = hd(x(list(B,:)),q);
end
xboot = sort(xboot);
lo = round(nboot*(alpha/2));
hi = nboot - lo;

% confidence intervals
CI(1) = xboot(lo+1);
CI(2) = xboot(hi);
