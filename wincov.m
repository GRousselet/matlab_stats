function [wcov,g]=wincov(x,y,percent)

% function [wv,g]=wincov(x,percent)
% returns the winsorized covariance of x & y as wcov
% returns number of winsorized observations at each tail of sample as g.
% x must be a vector.
% percent must be between 0 and 100; the function winsorizes the lower and
% upper extreme g values of x, where g=floor((percent/100)*length(x)).  If x is empty, then NaN is
% returned. 
% percent=20 by default
%
% GAR, University of Glasgow, Dec 2007

if nargin < 3;percent=20;end
if nargin < 2
    error('winvar:TooFewInputs', 'winvar requires two input arguments.');
elseif percent >= 100 || percent < 0
    error('winvar:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

if percent > 0 && percent < 1
    percent = percent * 100;
end
    
% make sure that x is a vector
sz = size(x);

dim=length(x);

% dim = find(sz == 1, 1);
if isempty(dim)    
  error('winvar:NonVectorData', 'winvar requires x to be a vector, not a matrix.');
end

[wx,g]=winsample(x,percent);
[wy,g]=winsample(y,percent);
tmp=cov(wx,wy);
wcov=tmp(1,2);

return
