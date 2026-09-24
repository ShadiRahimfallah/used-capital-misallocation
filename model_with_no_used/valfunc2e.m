function [v, vee, vwe] = valfunc2e(c, fspace, s)
global P Py fspacew cw sminw smaxw exit_prob VALFUNC2E_LOOP ...
    VALFUNC2E_CACHE_OFF

ns = size(s, 1);
iz = s(:,2);
iy = s(:,3);

Nz = size(P,  2);
Ny = size(Py, 2);

have_W = ~isempty(cw) && ~isempty(fspacew);

a_toW = min(max(s(:,1), sminw(1)), smaxw(1));

fast_ok = false;
try
    nE = fspace.n;
    pE = fspace.parms;
    fast_ok = numel(nE) == 3 ...
        && strcmp(fspace.bastype{2},'spli') && strcmp(fspace.bastype{3},'spli') ...
        && pE{2}{3} == 1 && pE{3}{3} == 1 ...
        && isequal(pE{2}{1}(:),(1:nE(2))') && isequal(pE{3}{1}(:),(1:nE(3))') ...
        && nE(2) == Nz && nE(3) == Ny;
    if fast_ok && have_W
        nW = fspacew.n;
        pW = fspacew.parms;
        fast_ok = numel(nW) == 3 ...
            && strcmp(fspacew.bastype{2},'spli') && strcmp(fspacew.bastype{3},'spli') ...
            && pW{2}{3} == 1 && pW{3}{3} == 1 ...
            && isequal(pW{2}{1}(:),(1:nW(2))') && isequal(pW{3}{1}(:),(1:nW(3))') ...
            && nW(2) == Nz && nW(3) == Ny;
    end
catch
    fast_ok = false;
end

force_loop = ~isempty(VALFUNC2E_LOOP) && logical(VALFUNC2E_LOOP);

if force_loop || ~fast_ok
    vee       = zeros(ns, 1);
    vwe       = zeros(ns, 1);
    v_survive = zeros(ns, 1);

    for kz = 1:Nz
        wz = P(iz, kz);
        for ky = 1:Ny
            wy  = Py(iy, ky);
            wzy = wz .* wy;

            ge = [s(:,1), kz*ones(ns,1), ky*ones(ns,1)];
            gw = [a_toW,  kz*ones(ns,1), ky*ones(ns,1)];

            ve_kzky = funeval_fast(c(:,1), fspace, ge(:,1), ge(:,2), ge(:,3));
            if have_W
                vw_kzky = funeval_fast(cw(:,1), fspacew, gw(:,1), gw(:,2), gw(:,3));
            else
                vw_kzky = -5 * ones(ns, 1);
            end

            vee       = vee + wzy .* ve_kzky;
            vwe       = vwe + wzy .* vw_kzky;
            v_survive = v_survive + wzy .* max(ve_kzky, vw_kzky);
        end
    end

    if have_W
        vdeath = vwe;
    else
        vdeath = -5 * ones(ns, 1);
    end

    ep = 0;
    if ~isempty(exit_prob), ep = exit_prob; end
    v = (1 - ep) .* v_survive + ep .* vdeath;
    v = v(:);
    return
end

persistent K_s1 K_atoW K_P K_Py pBaE pBaW pWgt K_cw pVW pvwe
cache_on = isempty(VALFUNC2E_CACHE_OFF) || ~VALFUNC2E_CACHE_OFF;

NaE = nE(1);

if ~(cache_on && isequal(K_s1,s(:,1)) && isequal(K_atoW,a_toW) ...
        && isequal(K_P,P) && isequal(K_Py,Py))
    pBaE   = splibas(pE{1}{1}, pE{1}{2}, pE{1}{3}, s(:,1));
    pBaW   = [];
    Pz     = P(iz,:);
    Pyy    = Py(iy,:);
    pWgt   = repmat(Pz,1,Ny) .* kron(Pyy, ones(1,Nz));
    K_s1   = s(:,1);
    K_atoW = a_toW;
    K_P    = P;
    K_Py   = Py;
    K_cw   = [];
end

CmE = reshape(full(c(:,1)), NaE, Nz*Ny);
VE  = pBaE * CmE;

if have_W
    if isempty(pBaW)
        pBaW = splibas(pW{1}{1}, pW{1}{2}, pW{1}{3}, a_toW);
    end
    if ~(cache_on && isequal(K_cw, cw(:,1)))
        pVW  = pBaW * reshape(full(cw(:,1)), nW(1), Nz*Ny);
        pvwe = sum(pWgt .* pVW, 2);
        K_cw = cw(:,1);
    end
    VW  = pVW;
    vwe = pvwe;
else
    VW  = -5 * ones(ns, Nz*Ny);
    vwe = sum(pWgt .* VW, 2);
end

vee       = sum(pWgt .* VE, 2);
v_survive = sum(pWgt .* max(VE, VW), 2);
vdeath    = vwe;

ep = 0;
if ~isempty(exit_prob), ep = exit_prob; end
v = (1 - ep) .* v_survive + ep .* vdeath;
v = v(:);
end
