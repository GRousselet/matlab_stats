function out = madr(x)
% out = madr(x)
% Computes the median absolute deviation to the median, 
% adjusting by a factor for asymptotically normal consistency.
% This is the formula used in the R version of mad, 
% which is different from the Matlab version.

% Copyright (C) 2012 Guillaume Rousselet - University of Glasgow

out = 1.4826.*median(abs(x-median(x)));