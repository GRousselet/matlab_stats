function d = cliffdelta(x,y)
% d = cliffdelta(x,y)
% Adaptation of Rand Wilcox's cid R function,
% from Rallfun-v19
% http://dornsife.usc.edu/labs/rwilcox/software/
% This version only returns delta.
% The full function is cid.
%
% INPUTS:
%         x & y are two vectors
%
% OUTPUTS:
%         d statistic P(X>Y)-P(X<Y)
%
% See also CID

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

x=x(:);
y=y(:);

yy=repmat(y,[1 length(x)])';
xx=repmat(x,[1 length(y)]);

m=xx-yy;
m=sign(m);
d=mean(m(:));

