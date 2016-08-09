function [mapd, apd] = mapd_eeg(x,y)
% [mpad, apd] = mapd_eeg(x,y)
% Computes all the pairwise differences between elements of x & y.
% Also estimates the median of these pairwise differences using
% the Harrell-Davis estimate of the 50th percentile.
% This statistics is a useful measure of effect size;
% it provides information about the typical difference between any two
% observations from two groups.
%
% INPUTS:
%         x & y are two matrices with dimensions time points x elements
%         elements could be for instance participants or trials
%
% OUTPUTS:
%         mapd is a vector of the medians of all the pairwise differences
%              dimensions = time points x 1
%         apd is a matrix time points x all the pairwise differences
%
% Adaptation of Rand Wilcox's wmwloc R function from Rallfun-v26
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, L2DCI, CID, CLIFFDELTA, MAPD

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow
% First version 2016-08-09
% to do: vectorise time dimension

[Nf,Nx] = size(x);
Ny = size(y,2);
mapd = zeros(Nf,1);
apd = zeros(Nf,Nx*Ny);

hx = x;
hy = y;

for F = 1:Nf

    x = hx(F,:);
    y = hy(F,:);

    % remove NaNs
    x = x(~isnan(x));
    y = y(~isnan(y));

    x = x(:);
    y = y(:);

    yy = repmat(y,[1 length(x)])';
    xx = repmat(x,[1 length(y)]);

    m = xx-yy; % up to this point, calculations are identical to Cliff's delta

    mapd(F) = hd(m(:));
    apd(F,:) = m;

end
