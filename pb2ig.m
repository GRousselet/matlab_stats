function [bootdiff,diff,CI,p,sig] = pb2ig(data1,data2,nboot,alpha,est,q,mu)

%pb2ig - percentile bootstrap test to compare 2 independent groups
% [bootdiff,diff,CI,p,sig] = pb2ig(data1,data2,nboot,alpha,est,q)
% 2 tailed bootstrap test based on the estimation of H1, the
% hypothesis of an experimental effect.
%
% INPUTS:
% -------
%
% DATA1 and DATA2 are 2 vectors:
% 
%           nboot (resamples, default 1000)
%           alpha (default 5%)
%           est (estimator, default 'median')
%           q (argument for estimator, default .5 for Harrell-Davis
%           estimator, 20% for trimmed mean
%           mu (optional null hypothesis, default 0)
%
% OUTPUTS:
% --------
%
% BOOTDIFF: bootstrapped distributions of differences
%
% LOCI/HICI: low and high boundaries of the confidence interval for
% the difference between 2 groups. Contrary to permtest and
% pbH0, the confidence interval calculated here is not under the
% null hypothesis. Therefore, instead of being centred on zero, it will
% follow the difference between the 2 experimental conditions compared.
% As a simple decisional rule, when the confidence interval under H1 does
% not include zero, the difference is considered significant.
%
% P value of the effect
%
% SIG: significativity of the experimental difference, binary output,
% 1=YES, 0=NO
%
% See also pbci bootse pb2dg hd tm

% Copyright (C) 2008, 2016 Guillaume Rousselet - University of Glasgow

% Guillaume A. Rousselet - University of Glasgow - January 2008
% changed randsample to randi - GAR - 25/05/2016

if nargin<3 || isempty(nboot)
  nboot=1000;
end
if nargin<4 || isempty(alpha)
    alpha=.05;
end
if nargin<5 || isempty(est)
est='median';
end

if nargin<6 || isempty(q)
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            q=.5;
        case {'tm'} % trimmed mean
            q=20;
        case {'trimmean'} % trimmed mean - Matlab function
            q=20; % later on multiplied by 2 to achieve 20% trimming of both tails
    end
end

n1 = length(data1);
n2 = length(data2);
lo = round(nboot.*alpha./2); % get CI boundary indices
hi = nboot - lo;
lo = lo+1;

for kk = 1:nboot % bootstrap with replacement loop
    switch lower(est)
    case {'hd'} % Harrell-Davis estimator
            eval(['bootdiff(kk)=',est,'(data1(randi(n1,n1,1)),q) -',est,'(data2(randi(n2,n2,1)),q);'])
    case {'tm'} % trimmed mean
            eval(['bootdiff(kk)=',est,'(data1(randi(n1,n1,1)),q) -',est,'(data2(randi(n2,n2,1)),q);'])
    case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootdiff(kk)=',est,'(data1(randi(n1,n1,1)),q.*2) -',est,'(data2(randi(n2,n2,1)),q.*2);'])
    otherwise
            eval(['bootdiff(kk)=',est,'(data1(randi(n1,n1,1))) -',est,'(data2(randi(n2,n2,1)));'])
    end    
end

diffsort = sort(bootdiff); % sort in ascending order
CI(1) = diffsort(lo);
CI(2) = diffsort(hi);

if nargout > 1
    if nargin < 7 || isempty(mu)
        mu=0;
    end
    p = sum(bootdiff<mu)./nboot+sum(bootdiff==mu)./(2.*nboot);
    p = 2.*min(p,1-p);
    sig = p<=alpha;
    
    switch lower(est)
        case {'hd'} % Harrell-Davis estimator
            eval(['diff=',est,'(data1,q) -',est,'(data2,q);'])
        case {'tm'} % trimmed mean
            eval(['diff=',est,'(data1,q) -',est,'(data2,q);'])
        case {'trimmean'} % trimmed mean - Matlab function
            eval(['diff=',est,'(data1,q.*2) -',est,'(data2,q.*2);'])
        otherwise
            eval(['diff=',est,'(data1) -',est,'(data2);'])
    end
    
end
 
      

       