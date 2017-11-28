function [Ty,CI,diff,se,df,p]=yuend(a,b,percent,alpha)

% function [Ty,CI,diff,se,df,p]=yuend(a,b,percent,alpha)
%
% Computes Ty (Yuen's T statistic) for two dependent groups
% Ty=(tma-tmb) / sqrt(da+db-2dab), where tma & tmb are trimmed means of a & b,
% da & db are yuen's estimate of the standard errors of tma & tmb.
% Ty is distributed approximately as Student's t with estimated degrees of freedom, df.
% The p-value (p) is the probability of obtaining a t-value whose absolute value
% is greater than Ty if the null hypothesis (i.e., tma-tmb = mu) is true.
% In other words, p is the p-value for a two-tailed test of H0: tma-tmb=0;
% Data arrays a & b must be vectors; percent must be a number between 0 & 100.
% se is the estimated standard error of the difference between the sample trimmed means. 
%
%   Default values:
%   mu = 0; 
%   percent = 20;
%   alpha = 0.05; necessary to compute the CI
%
% See Wilcox (2005), Introduction to Robust Estimation and Hypothesis
% Testing (2nd Edition), page 188-191 for a description of the Yuen
% procedure for dependent groups.
%
% GAR, University of Glasgow, Dec 2007
%
% See also YUEN 

if nargin<4;alpha=.05;end
if nargin<3;percent=20;end

if isempty(a) || isempty(b) 
    error('yuen:InvalidInput', 'data vectors cannot have length=0');
end

if (min(size(a))>1) || (min(size(b))>1)
    error('yuen:InvalidInput', 'yuen requires that the data are input as vectors.');
end

if (percent >= 100) || (percent < 0)
    error('yuen:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

[swa,ga]=winvar(a,percent); % winsorized variance of a & # items winsorized
[swb,gb]=winvar(b,percent); % winsorized variance of b & # items winsorized

% yuen's estimate of standard errors for a and b
na=length(a);
ha=na-2.*ga; % effective sample size after trimming
da=(na-1).*swa;

nb=length(b);
hb=nb-2.*gb;
db=(nb-1).*swb;

dab=(na-1).*wincov(a,b,percent);

% trimmed means
ma=tm(a,percent);
mb=tm(b,percent);

diff=ma-mb;

df= ha-1;
se=sqrt( (da+db-2.*dab)./(ha.*(ha-1)) );

Ty=(ma-mb)./se;

p=2*(1-tcdf(abs(Ty),df)); % 2-tailed probability

t=tinv(1-alpha./2,df); % 1-alpha./2 quantile of Student's distribution with df degrees of freedom
 
CI(1)=diff-t.*se; 
CI(2)=diff+t.*se;

return

% Data from Wilcox p.190
% before=[190 210 300 240 280 170 280 250 240 220];
% after=[210 210 340 190 260 180 200 220 230 200];

