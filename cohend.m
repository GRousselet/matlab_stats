function cod = cohend(x,y)
% cod = cohend(x,y);
% Computes Cohen's d for two independent groups using pooled standard deviation.
% x & y are two vectors

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% remove NaNs & reformat
x=x(~isnan(x)); x=x(:); n1 = numel(x);
y=y(~isnan(y)); y=y(:); n2 = numel(y);

diff = mean(x) - mean(y);
s1 = var(x,0);
s2 = var(y,0);
psd = sqrt( ((n1-1)*s1+(n2-1)*s2) / (n1+n2-2) ); % pooled standard deviation
cod = diff ./ psd;
