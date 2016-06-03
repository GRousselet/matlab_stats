function bse = bootse(x,nboot,est,q)

% bse = bootse(x,nboot,est,q)
% bootstrap estimate of the standard error of any estimator
%
% x = vector of observations
% nboot = number of bootstrap samples, e.g. 1000
% est = estimator, e.g. 'median', default 'hd'
% q = quantile for hd estimator
%
% See also pbci pb2ig pb2dg hd

% Copyright (C) 2007, 2013 Guillaume Rousselet - University of Glasgow

% GAR, University of Glasgow, Dec 2007
% GAR, April 2013: replace randsample with randi
% GAR, 2016-06-01: removed list option, edited help

if nargin<2;nboot=1000;est='hd';q=.5;end

x = x(:);
n = numel(x); % number of observations
boot = zeros(1,nboot);
list = randi(n,nboot,n); % bootstrap samples

switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(x(list(kk,:)),q);'])
        end
    otherwise
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(x(list(kk,:)));'])
        end
end

% We estimate the standard error of the test statistic by the 
% standard deviation of the bootstrap replications
% Efron & Tibshinari 1993, chapter 6
% Wilcox 2005, p.44-45
bse = std(boot,0); % normalize by (n-1)
