global W r Rnew qu pu zeta delta kappa phi gamma alpha eta theta F lambda_u ...
    egrid ygrid k Ny ce cw fspacee fspacew sminw smaxw smine smaxe exit_prob ...
    productivity_weight INV_ME INV_MW INV_AG INV_EXIT_RATE MKT_DEMAND_SUPPLY ...
    LAST_MOMENTS

if isempty(productivity_weight), productivity_weight = 1; end

compute_invariant_fast;

ag  = INV_AG(:);
S   = gridmake(ag, (1:k)', (1:Ny)');
a_s = S(:,1);
z_s = S(:,2);
y_s = S(:,3);

ne  = INV_ME(:);   nw = INV_MW(:);
tot = sum(ne) + sum(nw);
if tot > 0, ne = ne/tot;  nw = nw/tot; end

Eeff = productivity_weight .* exp(egrid(z_s));
[Kn_e, Ku_e, Kagg_e, L_e, Y_e, D_e, mu_e, bind_e0] = ...
    compute_capital_allocation(a_s, Eeff, Rnew, qu, phi, gamma, ...
    alpha, eta, W, F, pu, theta, 'vectorized', lambda_u);
Debt_e = max(0, Kn_e + pu .* Ku_e - a_s);

wage_w = W .* ygrid(y_s);
a_toW  = min(max(a_s, sminw(1)), smaxw(1));
VE_at  = funeval_fast(ce(:,1), fspacee, a_s,   z_s, y_s);
VW_at  = funeval_fast(cw(:,1), fspacew, a_toW, z_s, y_s);
p_exit = double(VW_at > VE_at);

ent_share_sim  = sum(ne);
L_demand_sim   = ne' * L_e;
L_supply_sim   = nw' * ygrid(y_s);
L_supply_heads = sum(nw);

Lbar_worker    = L_supply_sim / max(L_supply_heads, 1e-20);
empl_share_sim = sum(ne(L_e(:) >= Lbar_worker));

firm_size_eff   = L_demand_sim / max(ent_share_sim, 1e-20);
firm_size_heads = firm_size_eff / max(Lbar_worker, 1e-20);

Knew_demand_sim = ne' * Kn_e;
Dtot_sim        = ne' * a_s + nw' * a_s;
Y_total_sim     = ne' * Y_e;
Debt_total_sim  = ne' * Debt_e;
DtoY_sim        = Debt_total_sim / max(Y_total_sim, 1e-10);

Knew_supply_sim = Dtot_sim;
gap_L           = (L_demand_sim    - L_supply_sim)    / max(L_supply_sim,    1e-10);
gap_Kn          = (Knew_demand_sim - Knew_supply_sim) / max(Knew_supply_sim, 1e-10);

MKT_DEMAND_SUPPLY = struct('L_d',L_demand_sim, 'L_s',L_supply_sim, 'gap_L',gap_L, ...
    'Kn_d',Knew_demand_sim, 'Kn_s',Knew_supply_sim, 'gap_Kn',gap_Kn);

if exist('local_gini','file') == 2 && sum(nw) > 0
    gini_wage_model = local_gini(wage_w, nw);
else
    gini_wage_model = NaN;
end

gini_wage_trim = NaN;
if exist('local_gini','file') == 2 && sum(nw) > 0
    [ws_t, ord_t] = sort(wage_w(:));
    wm_t = nw(ord_t);
    st_  = sum(wm_t);
    if st_ > 0
        wm_t   = wm_t / st_;
        chi_   = cumsum(wm_t);
        clo_   = chi_ - wm_t;
        keepm_ = max(0, min(chi_, 0.99) - max(clo_, 0.01));
        if sum(keepm_) > 0
            gini_wage_trim = local_gini(ws_t, keepm_);
        end
    end
end

[vs, ord] = sort(L_e(:), 'descend');  ws = ne(ord);  Wt = sum(ws);
if Wt > 0
    cwm = cumsum(ws);  tgt = 0.10*Wt;  ix = find(cwm >= tgt, 1, 'first');
    if isempty(ix), ix = numel(ws); end
    frac              = min(max((tgt - (cwm(ix)-ws(ix))) / max(ws(ix),1e-20), 0), 1);
    top10_share_model = (sum(ws(1:ix-1).*vs(1:ix-1)) + frac*ws(ix)*vs(ix)) / max(sum(ws.*vs),1e-20);
else
    top10_share_model = NaN;
end

if ~isempty(INV_EXIT_RATE) && isfinite(INV_EXIT_RATE)
    firm_exit_model = INV_EXIT_RATE;
else
    ep = exit_prob;
    if isempty(ep), ep = 0;
    end
    if sum(ne) > 0
        firm_exit_model = ep + (1-ep) * (sum(ne .* p_exit) / sum(ne));
    else
        firm_exit_model = NaN;
    end
end

bind_e       = double(bind_e0);
constr_share = sum(ne .* bind_e) / max(sum(ne), 1e-20);
[~, ordC] = sort(L_e(:), 'descend');  wsC = ne(ordC);  bsC = bind_e(ordC);  muC = mu_e(ordC);
WtC = sum(wsC);
if WtC > 0
    ixC = find(cumsum(wsC) >= 0.10*WtC, 1, 'first');
    if isempty(ixC), ixC = numel(wsC); end
    constr_top10 = sum(wsC(1:ixC).*bsC(1:ixC)) / max(sum(wsC(1:ixC)), 1e-20);
    mu_top10     = sum(wsC(1:ixC).*muC(1:ixC)) / max(sum(wsC(1:ixC)), 1e-20);
else
    constr_top10 = NaN;
    mu_top10     = NaN;
end

a_top_    = max(a_s);
frac_top_ = sum(ne(a_s >= 0.999*a_top_)) + sum(nw(a_s >= 0.999*a_top_));
if frac_top_ > 0.005, sat_ = '  <-- SATURATING: raise amax'; else, sat_ = ''; end
fprintf('  [amax-check] top wealth node a=%.2f  mass@top=%.3f%%%s\n', a_top_, 100*frac_top_, sat_);

fprintf(['  [moments_from_invariant] ent=%.3f  empl=%.3f (L>=1 avg wkr)  gapL=%+.1f%%  gapKn=%+.1f%%  DtoY=%.3f\n' ...
    '      wGini=%.3f (trim %.3f)  top10emp=%.3f  exit=%.3f   [single capital: Ku=0]\n'], ...
    ent_share_sim, empl_share_sim, 100*gap_L, 100*gap_Kn, DtoY_sim, ...
    gini_wage_model, gini_wage_trim, top10_share_model, firm_exit_model);
fprintf('      constrained: binding share=%.3f  |  TOP-10%% firms by L: bind=%.3f  mean mu=%.3f\n', ...
    constr_share, constr_top10, mu_top10);
fprintf('      firm size: %.2f AVG WORKERS/firm  (%.3f efficiency units; Lbar/worker=%.3f)\n', ...
    firm_size_heads, firm_size_eff, Lbar_worker);

LAST_MOMENTS = struct('top10',top10_share_model, ...
    'ent',ent_share_sim,'empl',empl_share_sim,'exit',firm_exit_model,'wGini',gini_wage_model, ...
    'wGini_trim',gini_wage_trim, ...
    'DtoY',DtoY_sim,'constr',constr_share,'constr_top10',constr_top10, ...
    'firm_size_heads',firm_size_heads,'firm_size_eff',firm_size_eff,'Lbar_worker',Lbar_worker, ...
    'Y',Y_total_sim,'K',Knew_demand_sim,'KtoY',Knew_demand_sim/max(Y_total_sim,1e-10));
