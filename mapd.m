function [out, apd] = mapd(x,y)
% out = mapd(x,y)
% Estimates the median of the distribution of x-y.
% Returns the Harrell-Davis estimate of the median of all pairwise
% differences. This statistics is a useful measure of effect size;
% it provides information about the typical difference between any two
% observations from two groups.
%
% INPUTS:
%         x & y are two vectors
%
% OUTPUTS:
%         out is the median of all pairwise differences
%         apd is a vector containing all the pairwise differences
%
% Adaptation of Rand Wilcox's wmwloc R function from Rallfun-v26
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, L2DCI, CID, CLIFFDELTA

% Copyright (C) 2014, 2016 Guillaume Rousselet - University of Glasgow

% remove NaNs
x = x(~isnan(x));
y = y(~isnan(y));

x = x(:);
y = y(:);

yy = repmat(y,[1 length(x)])';
xx = repmat(x,[1 length(y)]);

apd = xx-yy; % up to this point, calculations are identical to Cliff's delta

out = hd(apd(:)); 
