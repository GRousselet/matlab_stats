function [phat,zhat] = disker(x,y,z)
% [phat,zhat] = disker(x,y,z)
% Estimate apparent effect size using the probability of correct 
% classification based on values in first group. 
% Uses adaptive kernel estimate with initial estimate based
% on expected frequency curve.
%
% x & y are two independent samples
% z is a vector of values classified as belonging or not to x
% If not specified, z = x.
% 
% phat = proportion of correctly classified observations
%   A "correct" classification is the event of deciding
%   that an observation from the first group did indeed 
%   come from the first group based on a kernel density 
%   estimate of the distributions.
%
% zhat = vector of 0s and 1s indicating  whether values 
%        in z would be classified as coming from group x.
% 
% Based on Rand Wilcox's disker R function
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also AKERD

% Copyright (C) 2012 Guillaume Rousselet - University of Glasgow

% remove NaNs & reformat
x=x(~isnan(x)); x=x(:); 
y=y(~isnan(y)); y=y(:); 

xsort = sort(x);

xhat = akerd(x,xsort,0);
yhat = akerd(y,xsort,0);%figure;hold on;plot(xsort,xhat,'k',xsort,yhat,'g')

% Compute apparent probability of a correct classification
phat = sum(xhat>yhat)./length(x);

if nargout==2
    
    if nargin<3 || isempty(z)
        z=x;
    end
    
    % Make decisions for the data in z,
    % set zhat=1 if decide it came from group 1.
    zxhat = akerd(x,z);
    zyhat = akerd(y,z);
    
    zhat = ones(length(z),1);
    zhat(zxhat<zyhat) = 0;
end
