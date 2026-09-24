function [v1, v2, x, v12, v21] = saveBelmaxw(c, fspace, s)
global beta fspacee ce smine smaxe P Py

ns = size(s, 1);

if size(c, 2) == 1 && length(c) == 2*ns
    c = [c(1:ns), c(ns+1:2*ns)];
end

solvew;

v1 = valfunc1w(c, fspace, s, x);
v2 = valfunc2w(c, fspace, s);

if nargout > 3
    v12 = beta * funbas(fspace, [x, s(:,2:end)]);

    Nz           = size(P,2);
    Ny           = size(Py,2);
    iz           = s(:,2);
    iy           = s(:,3);
    a            = s(:,1);
    a_after_clip = min(max(a, smine(1)), smaxe(1));
    feasible     = (a > smine(1));
    have_E       = ~isempty(ce) && ~isempty(fspacee);

    v21 = sparse(ns, ns);
    for kz = 1:Nz
        wz = P(iz, kz);
        for ky = 1:Ny
            wzy = wz .* Py(iy, ky);
            gw  = [a, kz*ones(ns,1), ky*ones(ns,1)];
            Phi = funbas(fspace, gw);
            VW  = Phi * c(:,1);
            if have_E
                ge = [a_after_clip, kz*ones(ns,1), ky*ones(ns,1)];
                VE = funeval(ce(:,1), fspacee, ge);
            else
                VE = -1e6 * ones(ns,1);
            end
            VE(~feasible) = -1e10;
            stay          = double(VW >= VE);
            v21           = v21 + bsxfun(@times, wzy.*stay, Phi);
        end
    end
end
end
