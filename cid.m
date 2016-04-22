function [qxly,q0,qxgy,d,ci] = cid(x,y,alpha)
% [qxly,q0,qxgy,d,ci] = cid(x,y,alpha)
% Adaptation of Wilcox's cid R function,
% from Rallfun-v19
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% Compute a confidence interval for delta using the method in
% Cliff, 1996, p. 140, eq 5.12
% Ordinal methods for behavioral data analysis. 
% Erlbaum, Mahwah, N.J.
%
% The null hypothesis is that for two independent group, P(X<Y)=P(X>Y).
% This function reports a 1-alpha confidence interval for
% P(X>Y)-P(X<Y)
%
% INPUTS:
%         x & y are two vectors
%         alpha is used to compute a (1-alpha) confidence interval
%
% OUTPUTS:
%         qxly: P(X<Y)
%         q0:   P(X=Y)
%         qxgy: P(X>Y)
%         d:    delta statistic P(X>Y)-P(X<Y)
%         ci:   confidence interval
%
% See also BINOMCI, CLIFFDELTA

% Copyright (C) 2012, 2016 Guillaume Rousselet - University of Glasgow
%
% ********************************************************************
% Example data set to compare to R output:
% x=[0 32 9 0 2 0 41 0 10 12 6 18 3 3 0 11 11 2 0 11];
% y=[1 2 3 2 5 6 6 1 8 0 3 0 2 10 32 12 2 0 3 0 5 8 3 2 4];
% xsort=sort(x);ysort=sort(y);
% xhat=akerd(x,xsort,0);yhat=akerd(y,ysort,0);
% figure;hold on;plot(x,zeros(size(x)),'kx',y,zeros(size(y)),'ro',xsort,xhat,'k',ysort,yhat,'r')
% x<-c(0,32,9,0,2,0,41,0,10,12,6,18,3,3,0,11,11,2,0,11)
% y<-c(1,2,3,2,5,6,6,1,8,0,3,0,2,10,32,12,2,0,3,0,5,8,3,2,4)

% remove NaNs
x=x(~isnan(x));
y=y(~isnan(y));

dx=length(x);
dy=length(y);

x=x(:);
y=y(:);

yy=repmat(y,[1 length(x)])';
xx=repmat(x,[1 length(y)]);

m=xx-yy;
msave=m;
m=sign(m);
d=mean(m(:));
phat=(1-d)./2;
flag=1;
if phat==0 || phat==1
    flag=0;
end
q0=sum(sum(msave==0))./numel(msave);
qxly=sum(sum(msave<0))./numel(msave);
qxgy=sum(sum(msave>0))./numel(msave);

if flag==1
    sigdih=sum(sum((m-d).^2))/(dx*dy-1);
    di=NaN(dx,1);
    for c=1:dx
        di(c)=sum(x(c)>y)/dy-sum(x(c)<y)./dy;
    end
    dh=NaN(dy,1);
    for c=1:dy
        dh(c)=sum(y(c)>x)/dx-sum(y(c)<x)/dx;
    end
    sdi=var(di);
    sdh=var(dh);
    sh=((dy-1)*sdi+(dx-1)*sdh+sigdih)/(dx*dy);
    zv=icdf('norm',alpha./2,0,1);
    cu=(d-d^3-zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh);
    cl=(d-d^3+zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh);
end

if flag==0
    nm=max(dx,dy);
    if phat==1
        ci=binomci(nm,nm,alpha);
    end
    if phat==0
        ci=binomci(0,nm,alpha);
    end
elseif flag==1
    ci=[cl cu];
end


