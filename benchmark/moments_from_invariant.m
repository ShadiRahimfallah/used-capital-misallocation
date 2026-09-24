global W r Rnew qu pu zeta delta kappa phi gamma alpha eta theta F lambda_u ...
    egrid ygrid k Ny ce cw fspacee fspacew sminw smaxw smine smaxe exit_prob ...
    productivity_weight INV_ME INV_MW INV_AG INV_EXIT_RATE MKT_DEMAND_SUPPLY ...
    DtoY KnKu_ratio beta_size_data wage_gini_target firm_exit_target ...
    top10_emp_share_target WGINI_USE_TRIM LAST_MOMENTS

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

Lbar_worker      = L_supply_sim / max(L_supply_heads, 1e-20);
empl_share_sim   = sum(ne(L_e(:) >= Lbar_worker));
Knew_demand_sim  = ne' * Kn_e;
Kused_demand_sim = ne' * Ku_e;
Dtot_sim         = ne' * a_s + nw' * a_s;
KnKu_sim         = Knew_demand_sim / max(Kused_demand_sim, 1e-10);
Y_total_sim      = ne' * Y_e;
Debt_total_sim   = ne' * Debt_e;
DtoY_sim         = Debt_total_sim / max(Y_total_sim, 1e-10);
Knew_supply_sim  = Dtot_sim / (1 + zeta*qu / (Rused*(delta+kappa)));
Kused_supply_sim = zeta * Knew_supply_sim / (delta + kappa);
gap_L            = (L_demand_sim     - L_supply_sim)     / max(L_supply_sim,     1e-10);
gap_Kn           = (Knew_demand_sim  - Knew_supply_sim)  / max(Knew_supply_sim,  1e-10);
gap_Ku           = (Kused_demand_sim - Kused_supply_sim) / max(Kused_supply_sim, 1e-10);

MKT_DEMAND_SUPPLY = struct('L_d',L_demand_sim, 'L_s',L_supply_sim, 'gap_L',gap_L, ...
    'Kn_d',Knew_demand_sim, 'Kn_s',Knew_supply_sim, 'gap_Kn',gap_Kn, ...
    'Ku_d',Kused_demand_sim,'Ku_s',Kused_supply_sim,'gap_Ku',gap_Ku);

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

us  = Ku_e ./ max(Kn_e + Ku_e, 1e-20);
okb = (ne > 0) & isfinite(us) & (a_s > 0);
if sum(okb) >= 5 && any(Ku_e(okb) > 1e-12)
    wb              = ne(okb)/sum(ne(okb));
    xb              = log(a_s(okb));
    yb              = us(okb);
    xbm             = sum(wb.*xb);
    ybm             = sum(wb.*yb);
    beta_size_model = sum(wb.*(xb-xbm).*(yb-ybm)) / max(sum(wb.*(xb-xbm).^2), 1e-20);
else
    beta_size_model = 0;
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

fprintf('\n');
fprintf('[moments_from_invariant]\n');

fprintf('\n');
fprintf('MARKET CLEARING\n');
fprintf('\n');
fprintf('gap Labor                     = %+8.1f%% \n', 100*gap_L);
fprintf('gap Knew                      = %+8.1f%% \n', 100*gap_Kn);

fprintf('\n');
fprintf('MOMENTS -- TARGETED\n');
fprintf('\n');
disp('                                1st col: model,    2nd col: data');
fprintf('\n');
fprintf('DtoY                          = %9.4f  %9.4f \n', DtoY_sim,          DtoY);
fprintf('KnKu                          = %9.4f  %9.4f \n', KnKu_sim,          KnKu_ratio);

if ~isempty(WGINI_USE_TRIM) && WGINI_USE_TRIM
    fprintf('wage Gini (trimmed, SCORED)   = %9.4f  %9.4f \n', gini_wage_trim,  wage_gini_target);
else
    fprintf('wage Gini (untrimmed, SCORED) = %9.4f  %9.4f \n', gini_wage_model, wage_gini_target);
end
fprintf('top10 employment share        = %9.4f  %9.4f \n', top10_share_model, top10_emp_share_target);
fprintf('exit rate                     = %9.4f  %9.4f \n', firm_exit_model,   firm_exit_target);
fprintf('beta_size                     = %9.4f  %9.4f \n', beta_size_model,   beta_size_data);

fprintf('\n');
fprintf('MOMENTS -- UNTARGETED\n');
fprintf('\n');
fprintf('entrepreneur share            = %9.4f \n', ent_share_sim);
fprintf('employment share (L>=1 wkr)   = %9.4f \n', empl_share_sim);
fprintf('wage Gini (untrimmed)         = %9.4f \n', gini_wage_model);

fprintf('\n');
fprintf('COLLATERAL CONSTRAINT\n');
fprintf('\n');
fprintf('binding share (all firms)     = %9.4f \n', constr_share);
fprintf('binding share (top-10%% by L)  = %9.4f \n', constr_top10);
fprintf('mean mu (top-10%% by L)        = %9.4f \n', mu_top10);

fprintf('\n');
fprintf('GRID DIAGNOSTIC\n');
fprintf('\n');
fprintf('top wealth node a             = %9.2f \n', a_top_);
fprintf('mass at top node              = %8.3f%%   %s \n', 100*frac_top_, sat_);
fprintf('\n');

LAST_MOMENTS = struct('top10',top10_share_model, ...
    'ent',ent_share_sim,'empl',empl_share_sim,'exit',firm_exit_model,'wGini',gini_wage_model, ...
    'wGini_trim',gini_wage_trim, 'beta_size',beta_size_model, ...
    'DtoY',DtoY_sim,'KnKu',KnKu_sim,'constr',constr_share,'constr_top10',constr_top10, ...
    'Y',Y_total_sim,'Kn',Knew_demand_sim,'Ku',Kused_demand_sim, ...
    'KntoY',Knew_demand_sim/max(Y_total_sim,1e-10), ...
    'KutoY',Kused_demand_sim/max(Y_total_sim,1e-10), ...
    'KtoY',(Knew_demand_sim+Kused_demand_sim)/max(Y_total_sim,1e-10));
