function [BEST,EST,CI,p] = pbci(x,Nb,alpha,est,q,mu)

% function [BEST,EST,CI,p] = pbci(x,nboot,alpha,est,q,mu)
%
% Computes a percentile bootstrap confidence interval
%
%   INPUTS:
%           x (vector)
%           Nb (resamples, default 1000)
%           alpha (default 5%)
%           est (estimator, default 'median')
%           q (argument for estimator, default .5 for Harrell-Davis
%           estimator, 20% for trimmed mean
%           mu (optional null hypothesis, necessary to get a p-value in
%           one-sample hypothesis testing procedure, default 0)
%
%   OUTPUTS:
%           BEST - bootstrap estimates
%           EST - estimate
%           CI - confidence interval
%           p - bootstrap p value
%
% See also bootse pb2ig pb2dg  

% Copyright (C) 2007, 2008, 2012 Guillaume Rousselet - University of Glasgow

% GAR - University of Glasgow - Dec 2007
% added EST output - Nov 2008
% added BEST output - Feb 2012

if nargin<2 || isempty(Nb)
    Nb=1000;
end
if nargin<3 || isempty(alpha)
    alpha=.05;
end
if nargin<4 || isempty(est)
    est='median';
end

if nargin<5 || isempty(q)
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            q=.5;
        case {'tm'} % trimmed mean
            q=20;
        case {'trimmean'} % trimmed mean - Matlab function
            q=20; % later on multiplied by 2 to achieve 20% trimming of both tails
    end
end

switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        eval(['EST=',est,'(x,q);'])
    case {'tm'} % trimmed mean
        eval(['EST=',est,'(x,q);'])
    case {'trimmean'} % trimmed mean - Matlab function
        eval(['EST=',est,'(x,q.*2);'])
    otherwise
        eval(['EST=',est,'(x);'])
end

n = length(x);
lo = round(Nb.*alpha./2);
hi = Nb - lo;
lo = lo+1;

for B = 1:Nb % bootstrap with replacement loop
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            eval(['bootx(B)=',est,'(x(randi(n,n,1)),q);'])
        case {'tm'} % trimmed mean
            eval(['bootx(B)=',est,'(x(randi(n,n,1)),q);'])
        case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootx(B)=',est,'(x(randi(n,n,1)),q.*2);'])
        otherwise
            eval(['bootx(B)=',est,'(x(randi(n,n,1)));'])
    end
end

diffsort = sort(bootx); % sort in ascending order
CI(1) = diffsort(lo);
CI(2) = diffsort(hi);
BEST = bootx;

if nargout > 2
    if nargin < 6 || isempty(mu)
        mu=0;
    end
    p = mean(bootx > mu) + mean(bootx == mu).*0.5;
    p = 2.*min(p,1-p);
end
