function [v, vew, vww] = valfunc2w(c, fspace, s)
global P Py fspacee ce smine smaxe VALFUNC2W_LOOP

ns = size(s, 1);
iz = s(:,2);
iy = s(:,3);

Nz = size(P,  2);
Ny = size(Py, 2);

have_E = ~isempty(ce) && ~isempty(fspacee);

a_after_entry = s(:,1);
feasible      = (a_after_entry > smine(1));
a_after_clip  = min(max(a_after_entry, smine(1)), smaxe(1));

fast_ok = false;
try
    nW = fspace.n;
    pW = fspace.parms;
    fast_ok = numel(nW) == 3 ...
        && strcmp(fspace.bastype{2},'spli') && strcmp(fspace.bastype{3},'spli') ...
        && pW{2}{3} == 1 && pW{3}{3} == 1 ...
        && isequal(pW{2}{1}(:),(1:nW(2))') && isequal(pW{3}{1}(:),(1:nW(3))') ...
        && nW(2) == Nz && nW(3) == Ny;
    if fast_ok && have_E
        nE = fspacee.n;
        pE = fspacee.parms;
        fast_ok = numel(nE) == 3 ...
            && strcmp(fspacee.bastype{2},'spli') && strcmp(fspacee.bastype{3},'spli') ...
            && pE{2}{3} == 1 && pE{3}{3} == 1 ...
            && isequal(pE{2}{1}(:),(1:nE(2))') && isequal(pE{3}{1}(:),(1:nE(3))') ...
            && nE(2) == Nz && nE(3) == Ny;
    end
catch
    fast_ok = false;
end

force_loop = ~isempty(VALFUNC2W_LOOP) && logical(VALFUNC2W_LOOP);

if force_loop || ~fast_ok
    vew = zeros(ns, 1);
    vww = zeros(ns, 1);
    v   = zeros(ns, 1);

    for kz = 1:Nz
        wz = P(iz, kz);
        for ky = 1:Ny
            wy  = Py(iy, ky);
            wzy = wz .* wy;

            gw = [s(:,1),       kz*ones(ns,1), ky*ones(ns,1)];
            ge = [a_after_clip, kz*ones(ns,1), ky*ones(ns,1)];

            vw_kzky = funeval_fast(c(:,1), fspace, gw(:,1), gw(:,2), gw(:,3));
            if have_E
                ve_kzky = funeval_fast(ce(:,1), fspacee, ge(:,1), ge(:,2), ge(:,3));
            else
                ve_kzky = -1e6 * ones(ns, 1);
            end

            ve_kzky(~feasible) = -1e10;

            vww = vww + wzy .* vw_kzky;
            vew = vew + wzy .* ve_kzky;
            v   = v   + wzy .* max(vw_kzky, ve_kzky);
        end
    end
    return
end

NaW = nW(1);
BaW = splibas(pW{1}{1}, pW{1}{2}, pW{1}{3}, s(:,1));
CmW = reshape(full(c(:,1)), NaW, Nz*Ny);
VW  = BaW * CmW;

if have_E
    NaE = nE(1);
    BaE = splibas(pE{1}{1}, pE{1}{2}, pE{1}{3}, a_after_clip);
    CmE = reshape(full(ce(:,1)), NaE, Nz*Ny);
    VE  = BaE * CmE;
else
    VE  = -1e6 * ones(ns, Nz*Ny);
end

VE(~feasible, :) = -1e10;

Pz  = P(iz,  :);
Pyy = Py(iy, :);
Wgt = repmat(Pz, 1, Ny) .* kron(Pyy, ones(1, Nz));

vww = sum(Wgt .* VW, 2);
vew = sum(Wgt .* VE, 2);
v   = sum(Wgt .* max(VW, VE), 2);
end
