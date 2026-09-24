function g = local_gini(x, w)
x = x(:);
if nargin < 2 || isempty(w)
    w = ones(size(x));
else
    w = w(:);
end
keep = isfinite(x) & isfinite(w) & w > 0 & x >= 0;
x    = x(keep);
w    = w(keep);
n    = numel(x);
if n == 0 || sum(w) == 0 || sum(x.*w) == 0
    g = NaN;
    return
end
[x, idx] = sort(x);
w     = w(idx);
cumw  = cumsum(w);
cumxw = cumsum(x .* w);
g     = 1 - sum( (cumxw + [0; cumxw(1:end-1)]) .* w ) / (cumxw(end) * sum(w));
end
