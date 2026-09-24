function v = valfunc1e(c, fspace, s, x)
global beta

v = menufun('fe', s, x) + beta * funeval_fast(c(:,2), fspace, x, s(:,2), s(:,3));
end
