function res = hd3d(x, q, nboot, alpha)
% res = HD3D(x,q,nboot, alpha)
% HD3D computes the 95% percentile bootstrap confidence interval
% for the qth quantile of a distribution using the Harrell-Davis estimator.
% The function uses a percentile bootstrap technique,
% which does not rely on an estimation of the standard error of hd.
% Unlike HD(), HD3D accepts a vector, 2D or 3D matrix as input.
%
% INPUTS:
% x Vector, 2D matrix or 3D matrix. HD is computed along the last
%   non-singleton dimension.
% q Quantile (default 0.5)
% nboot Number of bootstrap samples. If not specified confidence intervals
%       are not computed
% alpha Alpha level. Default to 0.05 if nboot is provided.
%
% OUTPUTS:
% res A structure with fields:
%   - xhd A vector or matrix of quantile estimates
%   - ci A vector or matrix of confidence intervals matching xhd
%
% EXAMPLES:
% res = hd3d(x) % compute HD(0.5) without confidence interval
% res = hd3d(x, 0.25) % specify we want the first quartile
% res = hd3d(x, 0.75, 2000) % estimate 3rd quartile & its 95% confidence interval
%                            using 2000 bootstrap samples
% res = hd3d(x, 1000, 0.10) % estimate 2nd quartile & its 90% confidence
%                             interval using 1000 bootstrap samples
%
% see:
% Wilcox, R.R. (2012)
% Introduction to robust estimation and hypothesis testing
% Academic Press
% p.126-132
%
% Adaptation of Rand Wilcox's qcipb R function,
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% See also HD, HDPBCI

% Copyright (C) 2017 Guillaume Rousselet - University of Glasgow
% GAR 2017-05-08 - first version

if nargin<2; q=0.5; end

x = squeeze(x);
if isvector(x)
    x = x(:)';
end
dim = ndims(x);
d = size(x);
n = d(dim);

% compute HD weights
m1 = (n+1).*q;
m2 = (n+1).*(1-q);
vec = 1:n;
w = betacdf(vec./n,m1,m2) - betacdf((vec-1)./n,m1,m2);

% compute quantiles
xsort = sort(x,dim);
if isvector(x)
    xhd = sum(w(:).*xsort(:));
else
    if dim == 2
        w2 = repmat(w,d(1),1);
        xhd = sum(w2 .* xsort, 2);
    elseif dim == 3
        w2 = permute(repmat(w,d(1),1,d(2)), [1 3 2]);
        xhd = sum(w2 .* xsort, 3);
    end
end
res.xhd = xhd;

% compute bootstrap confidence interval
if exist('nboot', 'var') && ~isempty(nboot)
    if ~exist('alpha', 'var') || isempty(alpha)
        alpha = 0.05;
    end
    lo = round(nboot*(alpha/2));
    hi = nboot - lo;
    lo = lo + 1;
    
    if isvector(x)
        xboot = sort(x(randi(n,nboot,n)),2);
        w2 = repmat(w,nboot,1);
        hdboot = sum(w2 .* xboot, 2);
        hdboot = sort(hdboot);
        ci(1) = hdboot(lo);
        ci(2) = hdboot(hi);
    else % handle matrix cases
        if dim == 2
            xboot = sort(reshape(x(:,randi(n,nboot,n)),[d(1), nboot, d(2)]), 3);
            w2 = permute(repmat(w,d(1),1,nboot), [1 3 2]);
            hdboot = sort(sum(w2 .* xboot, 3),2);
            ci(:,1) = hdboot(:,lo);
            ci(:,2) = hdboot(:,hi);
        elseif dim == 3
            bootsamples = randi(n,nboot,n);
            hdboot = zeros(d(1),d(2),nboot);
            w2 = permute(repmat(w,d(2),1,nboot), [1 3 2]);
            for E = 1:d(1)
            xboot = sort(reshape(squeeze(x(E,:,bootsamples)),[d(2), nboot, d(3)]), 3);
            hdboot(E,:,:) = sort(sum(w2 .* xboot, 3),2);
            end
            ci(:,:,1) = hdboot(:,:,lo);
            ci(:,:,2) = hdboot(:,:,hi);
        end
    end
    
    res.ci = ci;
    
end

%% check results and test function

% % 2d matrix
% rng(1)
% hx = randn(100,1);
% x = repmat(hx,1,10)'; % 10 x 100
% % 3d matrix
% x = zeros(2,10,100);
% x(1,:,:) = repmat(hx,1,10)';
% x(2,:,:) = repmat(hx,1,10)';
% 
% compare hd3d to hdpbci -------------------------
% nboot = 2000;
% alpha = 0.05;
% q = 0.5;
% 
% rng(1)
% x = randn(100,1);
% rng(1)
% [xhd,CI] = hdpbci(x,q,nboot,alpha);
% rng(1)
% res = hd3d(x,q,nboot,alpha);
% 
% rng(1)
% x = randn(10,100);
% hxhd = zeros(10,1);
% hCI = zeros(10,2);
% tic
% for d = 1:10
% rng(1)
% [xhd(d),hCI(d,:)] = hdpbci(x(d,:),q,nboot,alpha);
% end
% toc
% rng(1)
% tic
% res = hd3d(x,q,nboot,alpha);
% toc % ~50 times faster than the loop
% 
% compare to limo_harrell_davis -------------------------
% rng(1)
% x = randn(2,400,20); % electrodes x time points x participants
% rng(1)
% tic
% hdres = limo_harrell_davis(x,q,nboot);
% toc
% rng('default')
% rng(1)
% tic
% res = hd3d(x,q,nboot,alpha);
% toc % >400 times faster than limo_harrell_davis
% 
% 3d matrix -----------------------------------
% rng(1)
% hx = randn(100,1);
% x = zeros(2,10,100);
% x(1,:,:) = repmat(hx,1,10)';
% x(2,:,:) = repmat(hx,1,10)';
% hxhd = zeros(2,10);
% hCI = zeros(2,10,2);
% tic
% for E = 1:2
%     for F = 1:10
%         rng(1)
%         [xhd(E,F),hCI(E,F,:)] = hdpbci(squeeze(x(E,F,:)),q,nboot,alpha);
%     end
% end
% toc
% rng(1)
% tic
% res = hd3d(x,q,nboot,alpha);
% toc % >70 times faster than loop

