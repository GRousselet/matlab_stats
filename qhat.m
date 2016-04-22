function q = qhat(x,y,nboot)
% q = qhat(x,y,nboot)
% Estimates Q, a nonparametric measure of effect size for 
% two independent samples, as described in:
% Wilcox, R.R. & Muska, J. (2010) 
% Measuring effect size: A non-parametric analogue of omega(2). 
% The British journal of mathematical and statistical psychology, 52, 93-110.
%
% Uses an adaptive kernel estimate with an initial estimate based
% on the expected frequency curve.
% Bias is corrected using the .632 method of estimating prediction error
% (see Efron and Tibshirani, 1993, pp. 252--254).
%
% x & y are two vectors of observations
% nboot is the number of bootstrap samples - default 100
% 
% Based on Rand Wilcox's qhat R function
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also DISKER, AKERD

% Copyright (C) 2012, 2016 Guillaume Rousselet - University of Glasgow

x=x(~isnan(x));x=x(:);Nx=length(x);
y=y(~isnan(y));y=y(:);Ny=length(y);

if nargin<3 || isempty(nboot)
    nboot=100;
end

if nargin<4
    datax=randi(Nx,nboot,Nx);
    % datax is an nboot by n matrix containing subscripts for bootstrap sample
    % associated with first group.
    datay=randi(Ny,nboot,Ny);
    % datay is an nboot by m matrix containing subscripts for bootstrap sample
    % associated with second group.
end

%   fprintf('Taking bootstrap samples. Please wait.\n')
  
% Make bidx & bidy -------------------------------------------------------
% bidx is a n by nboot matrix. If the jth bootstrap sample from
% 1, ..., n contains the value i, bid[i,j]=0; otherwise bid[i,j]=1
bidx=zeros(Nx,nboot);
bidy=zeros(Ny,nboot);
for B = 1:nboot
    [bidx(:,B),~] = ismember(1:Nx,datax(B,:));
    [bidy(:,B),~] = ismember(1:Ny,datay(B,:));
end
bidx=bidx==0;
bidy=bidy==0;

temp3=zeros(nboot,Nx);
temp5=zeros(nboot,Ny);

for B=1:nboot
    [~,temp3(B,:)] = disker(x(datax(B,:)),y(datay(B,:)),x);
    % temp3 contains vector of 0s and 1s, 1 if x[i]
    % is classified as coming from group 1.
    [~,temp5(B,:)] = disker(y(datay(B,:)),x(datax(B,:)),y);
end

temp4=temp3.*bidx';
temp4=sum(temp4,1)./sum(bidx,2)';
temp6=temp5.*bidy';
temp6=sum(temp6,1)./sum(bidy,2)';
ep0x=nanmean(temp4); % epsilon hat_x
aperrorx=disker(x,y); % apparent error
regprex=.368.*aperrorx+.632.*ep0x;
ep0y=nanmean(temp6);
aperrory=disker(y,x); % apparent error
regprey=.368.*aperrory+.632.*ep0y;
q = (Nx.*regprex+Ny.*regprey)./(Nx+Ny);
