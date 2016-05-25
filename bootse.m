function bse = bootse(x,nboot,list,est,q)

% bse = bootse(x,nboot,list,est,q)
% bootstrap estimate of the standard error of any estimator
% 
% See also pbci pb2ig pb2dg 

% Copyright (C) 2007, 2013 Guillaume Rousselet - University of Glasgow

% GAR, University of Glasgow, Dec 2007
% GAR, April 2013: replace randsample with randi

if nargin<2;nboot=1000;est='hd';list=[];end
if nargin<5;q=.5;end

n=length(x); % number of subjects
boot = zeros(1,nboot);
if isempty(list)
   list=randi(n,nboot,n); 
end

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
