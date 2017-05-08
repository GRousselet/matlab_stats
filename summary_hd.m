function out = summary_hd(x, nboot, alpha)
% out = SUMMARY_HD(x, nboot, alpha)
% Compute a summary of observations stored in vector x:
% min, max and quartiles estimated using
% the Harrell-Davis estimator. Return results in a structure,
% as well as a table in the command window.
%
% INPUTS:
% x Vector of observations.
% nboot Number of bootstrap samples. If not specified confidence intervals
%       are not computed.
% alpha Alpha level. Default to 0.05 if nboot is provided.
%
% OUTPUTS:
% out A structure with fields:
%   - min Minimum observation(s)
%   - max Maximum observation(s)
%   - q Quantiles
%   - xq Quantile estimates
%   - qci Confidence intervals matching xhd
%
% EXAMPLES
% out = summary_hd(x) % basic call, does not return confidence intervals
% out = summary_hd(x, 2000) % specify number of bootstrap samples, alpha
%                             set to 0.05 by default
% out = summary_hd(x, 2000, 0.1) % specify 10% alpha level
%
% See also HD, HD3D

% Copyright (C) 2017 Guillaume Rousselet - University of Glasgow
% GAR 2017-05-08 - first version

if ~isvector(x)
    error('summary_hd only works with vectors...')
end

if ~exist('nboot', 'var')
    nboot = [];
end

if ~exist('alpha', 'var')
    alpha = [];
end

seq = [.25, .5, .75];
nq = numel(seq);

xmin = min(x);
xmax = max(x);
allq = zeros(nq,1);
if ~isempty(nboot)
    allci = zeros(nq,2);
end
for Q = 1:nq
    res = hd3d(x, seq(Q), nboot, alpha);
    allq(Q) = res.xhd;
    if ~isempty(nboot)
        allci(Q,:) = res.ci;
    end
end

out.min = xmin;
out.max = xmax;
out.q = seq;
out.xq = allq;
if ~isempty(nboot)
    out.qci = allci;
end

results = [xmax allq(3) allq(2) allq(1) xmin]';
if ~isempty(nboot)
    ci = [NaN NaN; allci(3,:); allci(2,:); allci(1,:); NaN NaN];
end

if ~isempty(nboot)
    T = table(results,ci, 'RowNames', {'maximum';'quartile 3';'quartile 2';'quartile 1';'minimum'});
else
    T = table(results, 'RowNames', {'maximum';'quartile 3';'quartile 2';'quartile 1';'minimum'});
end
T



