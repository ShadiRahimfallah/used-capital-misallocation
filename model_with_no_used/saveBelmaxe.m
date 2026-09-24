function [v1, v2, v3, x, v12, v21, v33] = saveBelmaxe(c, fspace, s)
global beta r fspacew cw exit_prob P Py sminw smaxw

ns = size(s, 1);

solvee;

v1 = valfunc1e(c, fspace, s, x);

[v2, vee, vwe] = valfunc2e(c, fspace, s);

v3 = zeros(size(v1));

if nargout > 4
    v12 = beta * funbas(fspace, [x, s(:,2:end)]);

    Nz    = size(P,2);
    Ny    = size(Py,2);
    iz    = s(:,2);
    iy    = s(:,3);
    a     = s(:,1);
    a_toW = min(max(a, sminw(1)), smaxw(1));
    ep    = 0;
    if ~isempty(exit_prob), ep = exit_prob;
    end
    have_W = ~isempty(cw) && ~isempty(fspacew);

    v21 = sparse(ns, ns);
    for kz = 1:Nz
        wz = P(iz, kz);
        for ky = 1:Ny
            wzy = wz .* Py(iy, ky);
            ge  = [a, kz*ones(ns,1), ky*ones(ns,1)];
            Phi = funbas(fspace, ge);
            VE  = Phi * c(:,1);
            if have_W
                gw = [a_toW, kz*ones(ns,1), ky*ones(ns,1)];
                VW = funeval(cw(:,1), fspacew, gw);
            else
                VW = -5 * ones(ns,1);
            end
            stay = double(VE >= VW);
            v21  = v21 + bsxfun(@times, (1-ep).*wzy.*stay, Phi);
        end
    end

    v33 = sparse(ns, ns);
end
end
