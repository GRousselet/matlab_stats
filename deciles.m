function thetaq = deciles(x)
% thetaq = deciles(x)
% Computes the Harrell-Davis estimates of the deciles.
% The vector x contains the data.
% Quantiles .1 .2 .3 .4 .5 .6 .7 .8 .9 are estimated
%
% Adaptation of Rand Wilcox's hd R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
% Original article:
% http://biomet.oxfordjournals.org/content/69/3/635.abstract

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

n = numel(x);
vec = 1:n;
y = sort(x);
thetaq = zeros(9,1);

for dec = 1:9
    
    q = dec/10;
    m1 = (n+1).*q;
    m2 = (n+1).*(1-q);
    w = betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
    thetaq(dec) = sum(w(:).*y(:));
    
end