function v = funeval_fast(cvec, fspace, ax, iz, iy)
persistent warned

ok = false;
try
    n = fspace.n;
    p = fspace.parms;
    ok = numel(n) == 3 ...
        && strcmp(fspace.bastype{2}, 'spli') && strcmp(fspace.bastype{3}, 'spli') ...
        && p{2}{3} == 1 && p{3}{3} == 1 ...
        && isequal(p{2}{1}(:), (1:n(2))') && isequal(p{3}{1}(:), (1:n(3))') ...
        && iscolumn(ax) && isequal(size(ax), size(iz), size(iy)) ...
        && all(iz == round(iz)) && all(iy == round(iy)) ...
        && all(iz >= 1 & iz <= n(2)) && all(iy >= 1 & iy <= n(3));
catch
    ok = false;
end
if ~ok
    if isempty(warned), warned = true;
        warning('funeval_fast: unexpected fspace/state layout -- falling back to funeval (slow).');
    end
    v = funeval(cvec, fspace, [ax, iz, iy]);
    return
end

Na   = n(1);
Nz   = n(2);
Ba   = splibas(p{1}{1}, p{1}{2}, p{1}{3}, ax);
Cm   = reshape(full(cvec), Na, []);
cols = iz + Nz * (iy - 1);
v    = full(sum(Ba .* (Cm(:, cols)).', 2));
end
