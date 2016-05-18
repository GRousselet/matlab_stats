function [est, ci, pval] = l2dci(x,y)
% [est, ci, pval] = l2dci(x,y)
% Computes a bootstrap confidence interval for a
% measure of location (hd) of the distribution of x-y.
% Adaptation of Rand Wilcox's l2dci R function from Rallfun-v26:
% http://dornsife.usc.edu/labs/rwilcox/software/
% 
% INPUTS:
%         x & y are two vectors
%
% OUTPUTS:
%         est  = Harrell-Davis estimate of the median of all pairwise differences
%               see mapd for details
%         ci   = confidence interval
%        pval  = p value
%
% See also HD, MAPD, CID

% Copyright (C) 2014, 2016 Guillaume Rousselet - University of Glasgow

alpha = .05;
Nb = 500;

% remove NaNs
x = x(~isnan(x));
y = y(~isnan(y));

nx = numel(x);
ny = numel(y);

est = mapd(x,y);

bvec = zeros(Nb,1);

for B = 1:Nb
    bvec(B) = mapd( x(randi(nx,nx,1)),y(randi(ny,ny,1)));
end

bvec = sort(bvec);
low = round((alpha/2)*Nb)+1;
up = Nb-low;
temp = sum(bvec<0)/Nb + sum(bvec==0)/(2*Nb);
pval = 2*(min(temp,1-temp));
ci = bvec([low up]);
% estimate of square standard error of the estimator used
% se = var(bvec)