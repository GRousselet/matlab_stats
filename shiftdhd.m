function [xd, yd, delta, deltaCI] = shiftdhd(x,y,nboot,plotit)
%[xd yd delta deltaCI] = shiftdhd(x,y,nboot,plotshift)
% SHIFTDHD computes the 95% simultaneous confidence intervals
% for the difference between the deciles of two dependent groups
% using the Harrell-Davis quantile estimator, in conjunction with
% percentile bootstrap estimation of the quantiles' standard errors.
% See Wilcox Robust hypothesis testing book 2005 p.184-187
%
% INPUTS:
% - x & y are vectors of the same length without missing values
% - nboot = number of bootstrap samples - default = 200
% - plotit = 1 to get a figure; 0 otherwise by default
%
% OUTPUTS:
% - xd & yd = vectors of quantiles
% - delta = vector of differences between quantiles (x-y)
% - deltaCI = matrix quantiles x low/high bounds of the confidence
% intervals of the quantile differences
%
% Tied values are not allowed.
% If tied values occur, or if you want to estimate quantiles other than the deciles,
% or if you want to use an alpha other than 0.05, use SHIFTDHD_PBCI instead.
%
% Adaptation of Rand Wilcox's shifthd R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, SHIFTDHD, SHIFTHD_PBCI

% Copyright (C) 2007, 2013, 2016 Guillaume Rousselet - University of Glasgow
% Original R code by Rand Wilcox
% GAR, University of Glasgow, Dec 2007
% GAR, April 2013: replace randsample with randi
% GAR, June 2013: added list option
% GAR, 7 Jun 2016: removed list option
% GAR, 17 Jun 2016: updated help for github release

if nargin < 3;nboot=200;end
if nargin < 4;plotit=0;end

if numel(x) ~= numel(y)
   error('shiftdhd: x and y must be the same length') 
end

n=length(x);
c=(37./n.^1.4)+2.75; % The constant c was determined so that the simultaneous
                     % probability coverage of all 9 differences is
                     % approximately 95% when sampling from normal
                     % distributions

% Get >>ONE<< set of bootstrap samples
% The same set is used for all nine quantiles being compared
list = randi(n,nboot,n);

xd = zeros(9,1);
yd = zeros(9,1);
delta = zeros(9,1);
deltaCI = zeros(9,2);
bootdelta = zeros(nboot,1);

for d=1:9
    q = d/10;
    xd(d) = hd(x,q);
    yd(d) = hd(y,q);
    delta(d) = xd(d) - yd(d);
    for b=1:nboot
        bootdelta(b) = hd(x(list(b,:)),q) - hd(y(list(b,:)),q);
    end
    delta_bse = std(bootdelta,0);
    deltaCI(d,1) = delta(d)-c.*delta_bse;
    deltaCI(d,2) = delta(d)+c.*delta_bse;
end

if plotit==1
    
    figure('Color','w','NumberTitle','off');hold on
    
    ext = 0.1*(max(xd)-min(xd));
    plot([min(xd)-ext max(xd)+ext],[0 0],'LineWidth',1,'Color',[.5 .5 .5]) % zero line
    for qi = 1:9
       plot([xd(qi) xd(qi)],[deltaCI(qi,1) deltaCI(qi,2)],'k','LineWidth',2) 
    end
    % mark median
    v = axis;plot([xd(5) xd(5)],[v(3) v(4)],'k:')
    plot(xd,delta,'ko-','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10,'LineWidth',1)
    set(gca,'FontSize',14,'XLim',[min(xd)-ext max(xd)+ext])
    box on
    xlabel('Group 1 deciles','FontSize',16)
    ylabel('Group 1 - group 2 deciles','FontSize',16)
    
end

%% TEST
% Data from Wilcox p.187
% time1=[0 32 9 0 2 0 41 0 0 0 6 18 3 3 0 11 11 2 0 11];
% time3=[0 25 10 11 2 0 17 0 3 6 16 9 1 4 0 14 7 5 11 14];
