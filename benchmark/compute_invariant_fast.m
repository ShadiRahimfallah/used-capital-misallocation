function compute_invariant_fast()
global cxe cxw ce cw fspacee fspacew P Py k Ny egrid sminw smaxw smine smaxe ...
    exit_prob productivity_weight LOCKOUT_EXOG_DEATH INV_NA INV_ME INV_MW ...
    INV_AG INV_USE_EIGS INV_CURV INV_EXIT_RATE COMPUTE_AGE INV_AGE_MASS ...
    INV_AGE_COH INV_AGE_KEEP

if isempty(LOCKOUT_EXOG_DEATH), LOCKOUT_EXOG_DEATH = true; end
if isempty(INV_USE_EIGS),       INV_USE_EIGS       = false; end
if ~isempty(INV_NA) && INV_NA > 10, Na = INV_NA; else, Na = 800; end
if isempty(INV_CURV) || INV_CURV <= 0, INV_CURV = 2.5; end

tol_inv = 1e-10;
max_inv = 20000;

amin = min(sminw(1), smine(1));
amax = max(smaxw(1), smaxe(1));

ag = amin + (amax - amin) * ( linspace(0, 1, Na)'.^INV_CURV );
ag = unique(ag(:));   Na = numel(ag);
Ns = Na*k*Ny;

S      = gridmake(ag, (1:k)', (1:Ny)');
a_s    = S(:,1);
z_s    = S(:,2);
y_s    = S(:,3);
zy_off = (z_s-1)*Na + (y_s-1)*Na*k;
ia_s   = mod((1:Ns)'-1, Na) + 1;

aw_eval = min(max(a_s, sminw(1)), smaxw(1));
ae_eval = min(max(a_s, smine(1)), smaxe(1));
VW      = funeval_fast(cw(:,1), fspacew, aw_eval, z_s, y_s);
VE      = funeval_fast(ce(:,1), fspacee, ae_eval, z_s, y_s);
xw      = min(max(funeval_fast(cxw, fspacew, aw_eval, z_s, y_s), sminw(1)), smaxw(1));
xe      = min(max(funeval_fast(cxe, fspacee, ae_eval, z_s, y_s), smine(1)), smaxe(1));

feas     = a_s > smine(1);
prefer_E = (VE >= VW);
canEnter = feas & prefer_E;

posE = min(max(interp1(ag, (1:Na)', xe, 'linear', 'extrap'), 1), Na);
loE  = min(floor(posE), Na-1);
wE   = posE - loE;
hiE  = loE + 1;
posW = min(max(interp1(ag, (1:Na)', xw, 'linear', 'extrap'), 1), Na);
loW  = min(floor(posW), Na-1);
wW   = posW - loW;
hiW  = loW + 1;
LloE = loE + zy_off;
LhiE = hiE + zy_off;
LloW = loW + zy_off;
LhiW = hiW + zy_off;

clipWval = min(max(ag, sminw(1)), smaxw(1));
jWclip   = min(max(round(interp1(ag, (1:Na)', clipWval, 'linear', 'extrap')), 1), Na);
aClipW   = jWclip(ia_s) + zy_off;

idx   = (1:Ns)';
SE    = sparse([LloE; LhiE], [idx; idx], [1-wE; wE], Ns, Ns);
SW    = sparse([LloW; LhiW], [idx; idx], [1-wW; wW], Ns, Ns);
Cclip = sparse(aClipW, idx, 1, Ns, Ns);

prefE   = double(prefer_E);
nprefE  = 1 - prefE;
cEnter  = double(canEnter);
ncEnter = 1 - cEnter;
p       = exit_prob;
lockout = LOCKOUT_EXOG_DEATH;

    function mu2 = Astep(mu)
        me_     = mu(1:Ns);
        mw_     = mu(Ns+1:end);
        me1     = inv_shock(me_, Na, k, Ny, P, Py);
        mw1     = inv_shock(mw_, Na, k, Ny, P, Py);
        survE   = (1-p)*me1;
        killed  = p*me1;
        Estay   = survE.*prefE;
        Eexit   = survE.*nprefE;
        Wenter  = mw1.*cEnter;
        Wstay   = mw1.*ncEnter;
        mLocked = Cclip*killed;
        mEexitW = Cclip*Eexit;
        if ~lockout
            Efinal = Estay + Wenter + mLocked.*cEnter;
            Wfinal = Wstay + mEexitW + mLocked.*ncEnter;
        else
            Efinal = Estay + Wenter;
            Wfinal = Wstay + mEexitW + mLocked;
        end
        mu2 = [SE*Efinal; SW*Wfinal];
    end

mu = [(0.25/Ns)*ones(Ns,1); (0.75/Ns)*ones(Ns,1)];

if INV_USE_EIGS
    opts.maxit = 1000;
    opts.tol = 1e-12;
    try
        [v, ~] = eigs(@(x) Astep(x), 2*Ns, 1, 'largestreal', opts);
    catch
        [v, ~] = eigs(@(x) Astep(x), 2*Ns, 1, 'lr', opts);
    end
    v = real(v);
    if sum(v) < 0, v = -v; end
    v(v < 0) = 0;
    mu       = v / sum(v);
    it       = 0;
    d        = 0;
else
    d = inf;
    for it = 1:max_inv
        mu_old = mu;
        mu     = Astep(mu);
        mu     = mu / sum(mu);
        d      = max(abs(mu - mu_old));
        if d < tol_inv, break; end
    end
end
me = mu(1:Ns);
mw = mu(Ns+1:end);

me1     = inv_shock(me, Na, k, Ny, P, Py);
mw1     = inv_shock(mw, Na, k, Ny, P, Py);
killed  = p*me1;
survE   = (1-p)*me1;
mLocked = Cclip*killed;
Estay   = survE.*prefE;
Eexit   = survE.*nprefE;
mEexitW = Cclip*Eexit;
Wenter  = mw1.*cEnter;
Wstay   = mw1.*ncEnter;
meP     = Estay + Wenter;
mwP     = Wstay + mEexitW + mLocked;

ent_incumbent = sum(me1);
if ent_incumbent > 0
    INV_EXIT_RATE = sum(killed + Eexit) / ent_incumbent;
else
    INV_EXIT_RATE = NaN;
end

INV_ME = meP;
INV_MW = mwP;
INV_AG = ag;
fprintf('  [compute_invariant_fast] converged (it=%d, d=%.1e, entShare=%.4f, Na=%d, eigs=%d)\n', ...
    it, d, sum(meP), Na, INV_USE_EIGS);

if ~isempty(COMPUTE_AGE) && COMPUTE_AGE
    AMAX  = 400;
    AKEEP = 40;
    if ~isempty(INV_AGE_KEEP), AKEEP = INV_AGE_KEEP; end
    age_mass = zeros(AMAX, 1);
    coh      = Wenter;
    INV_AGE_COH      = zeros(numel(coh), AKEEP);
    age_mass(1)      = sum(coh);
    INV_AGE_COH(:,1) = coh;
    for aa = 2:AMAX
        coh          = prefE .* ((1-p) * inv_shock(SE*coh, Na, k, Ny, P, Py));
        age_mass(aa) = sum(coh);
        if aa <= AKEEP, INV_AGE_COH(:,aa) = coh; end
        if age_mass(aa) < 1e-12 * max(age_mass(1), 1e-300), break; end
    end
    tot          = sum(age_mass);
    INV_AGE_MASS = age_mass;
    fprintf('  [firm-age] mean age = %.1f  frac aged 1-10 = %.4f  (cohort sum %.5f vs meP %.5f)\n', ...
        sum((1:AMAX)'.*age_mass)/max(tot,1e-300), ...
        sum(age_mass(1:min(10,AMAX)))/max(tot,1e-300), tot, sum(meP));
end
end
