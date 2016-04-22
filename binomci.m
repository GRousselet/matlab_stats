function bci = binomci(x,n,alpha)
% bci = binomci(x,n,alpha)
% Computes a 1-alpha confidence interval for p, the probability of
% success for a binomial distribution, using Pratt's method.
% 
% n is the number of trials
% x is the number of successes observed among n trials
% 
% Adaptation of Rand Wilcox's binomci R function,
% from Rallfun-v19
% http://dornsife.usc.edu/labs/rwilcox/software/

% Copyright (C) 2012 Guillaume Rousselet - University of Glasgow

if x~=n && x~=0
    z=icdf('norm',1-alpha/2,0,1);
    A=((x+1)./(n-x)).^2;
    B=81.*(x+1).*(n-x)-9.*n-8;
    C=(0-3).*z.*sqrt(9.*(x+1).*(n-x).*(9.*n+5-z.^2)+n+1);
    D=81.*(x+1).^2-9.*(x+1).*(2+z.^2)+1;
    E=1+A.*((B+C)/D).^3;
    upper=1./E;
    A=(x./(n-x-1)).^2;
    B=81.*x.*(n-x-1)-9.*n-8;
    C=3.*z.*sqrt(9.*x.*(n-x-1).*(9.*n+5-z.^2)+n+1);
    D=81.*x.^2-9.*x.*(2+z.^2)+1;
    E=1+A.*((B+C)./D).^3;
    lower=1./E;
end

if x==0
    lower=0;
    upper=1-alpha.^(1./n);
end
if x==1
    upper=1-(alpha./2).^(1./n);
    lower=1-(1-alpha./2).^(1./n);
end
if x==n-1
    lower=(alpha/2).^(1./n);
    upper=(1-alpha/2).^(1./n);
end
if x==n
    lower=alpha.^(1./n);
    upper=1;
end
bci=[lower upper];