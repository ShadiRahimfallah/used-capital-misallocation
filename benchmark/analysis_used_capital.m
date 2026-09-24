clear; clear global; clc;

global rho alpha eta beta delta W r smin smax fspacee fspacew fspace ce cw ...
    cxe cxw cse csw egrid P smine sminw smaxw theta xai F Lbar phi gamma se ...
    Rnew eps_model productivity_weight Perg k exit_prob smaxe qu zeta kappa ...
    pu use_calib_sim_size calib_mode DtoY KnKu_ratio ygrid Py Ny rho_y ...
    sigma_y Ey piy sigma_crra lambda_u sigma_perm eta_p omega ...
    fast_diagnostic beta_size_data wage_gini_target firm_exit_target ...
    top10_emp_share_target xsim_prev_for_amax amax_calib_mult ...
    hybrid_use_nested_clearing hybrid_inner_clearing_skip ...
    hybrid_nested_price_lb hybrid_nested_price_ub gaps_calibration_last ...
    LOCKOUT_EXOG_DEATH nested_clearing_quiet nested_clearing_max_iter ...
    nested_clearing_tol nested_clearing_tol_Knew calib_moment_weights ...
    USE_INVARIANT_SEED SIM_NF SIM_TINIT LAST_MOMENTS psi_u WGINI_USE_TRIM ...
    INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT INV_ME INV_MW INV_AG ...
    R_TARGET COMPUTE_AGE INV_AGE_MASS INV_AGE_COH INV_AGE_KEEP

sd = fileparts(mfilename('fullpath'));
cd(sd);
startup;

parameters_benchmark;

fast_diagnostic            = false;
LOCKOUT_EXOG_DEATH         = true;
calib_mode                 = true;
use_calib_sim_size         = true;
amax_calib_mult            = 1.9;
hybrid_use_nested_clearing = false;
hybrid_inner_clearing_skip = false;
USE_INVARIANT_SEED         = true;
SIM_NF                     = 12000;
SIM_TINIT                  = 2;
INV_NA                     = 3000;
INV_CURV                   = 2.5;
INV_USE_EIGS               = true;
MOMENTS_FROM_INVARIANT     = true;
nested_clearing_quiet      = true;
nested_clearing_tol        = 0.005;
nested_clearing_tol_Knew   = 0.005;
calib_moment_weights       = [30; 30; 30; 30; 30; 30];
WGINI_USE_TRIM             = false;
hybrid_nested_price_lb     = [0.35; 0.005; 0.04];
hybrid_nested_price_ub     = [2.00; 0.99*(1/beta-1); 0.26];

COMPUTE_AGE   = true;
INV_AGE_KEEP  = 40;

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic);
[ygrid, Py, Ey, piy]          = rouwenhorst_y(Ny, rho_y, sigma_y);
zeta = (delta + kappa) / KnKu_ratio;

theta_cal = theta;
W_cal     = W;
r_cal     = r;
Rnew_cal  = Rnew;

NFFILE = fullfile(sd, 'nofriction_result.mat');
if ~exist(NFFILE, 'file')
    error(['analysis_used_capital: %s not found. Run nofriction.m first -- this ' ...
        'analysis compares the benchmark against its own frictionless twin.'], NFFILE);
end
NF = load(NFFILE);
fn = fieldnames(NF);
NFR = NF.(fn{1});
if abs(NFR.base.theta - theta_cal) > 1e-6
    warning(['analysis_used_capital: nofriction_result.mat was built at theta=%.6f but ' ...
        'parameters_benchmark has theta=%.6f. The two economies are not the same point.'], ...
        NFR.base.theta, theta_cal);
end

ECON = { 'friction',     theta_cal, W_cal,             r_cal,             Rnew_cal
         'frictionless', 1e5,       NFR.frictionless.W, NFR.frictionless.r, NFR.frictionless.Rnew };
NE = size(ECON, 1);

fprintf('================================================================\n');
fprintf('  BENCHMARK vs ITS FRICTIONLESS TWIN -- life cycle and selection\n');
fprintf('  theta_cal=%.4f eta_p=%.4f omega=%.4f kappa=%.4f sigma_y=%.4f beta=%.4f\n', ...
    theta_cal, eta_p, omega, kappa, sigma_y, beta);
for i = 1:NE
    fprintf('  %-13s theta=%10.4g  W=%.4f r=%.4f Rnew=%.4f\n', ...
        ECON{i,1}, ECON{i,2}, ECON{i,3}, ECON{i,4}, ECON{i,5});
end
fprintf('  Prices are held at each economy''s own cleared equilibrium; nothing is re-cleared.\n');
fprintf('================================================================\n');

R = struct('name',{},'ushare',{},'mu',{},'MRPK',{},'me',{},'a',{},'z',{},'aq',{},'zt',{}, ...
    'age_used',{},'age_mu',{},'age_mass',{},'age_mrpk',{},'age_sdm',{},'tabA',{}, ...
    'sd_logMRPK',{},'ent',{},'W',{},'r',{},'Rnew',{},'qu',{},'pu',{},'theta',{});

for e = 1:NE
    theta = ECON{e,2};
    W     = ECON{e,3};
    r     = ECON{e,4};
    Rnew  = ECON{e,5};

    fprintf('\n\n----------------------------------------------------------------\n');
    fprintf('  SOLVING: %s   (theta=%.4g, W=%.4f, r=%.4f)\n', ECON{e,1}, theta, W, r);
    fprintf('----------------------------------------------------------------\n');

    x_full = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    for w = 1:3
        evaluate_economy(x_full);
    end
    evaluate_economy(x_full);

    [qu, pu, Rused] = derive_used_prices(r, delta, kappa, Rnew, zeta, psi_u);

    ag  = INV_AG(:);
    S   = gridmake(ag, (1:k)', (1:Ny)');
    a_s = S(:,1);
    z_s = S(:,2);
    me_raw = INV_ME(:);
    mw_raw = INV_MW(:);
    pop    = sum(me_raw) + sum(mw_raw);
    me_pop = me_raw / max(pop, 1e-300);
    mw_pop = mw_raw / max(pop, 1e-300);
    me     = me_raw / max(sum(me_raw), 1e-300);

    Eeff = productivity_weight .* exp(egrid(z_s));
    [Kn_e, Ku_e, Kagg_e, L_e, Y_e, ~, mu_e, bind_e] = compute_capital_allocation( ...
        a_s, Eeff, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, 'vectorized', lambda_u);

    MRPK      = (1 - alpha) * eta * Y_e ./ max(Kagg_e, 1e-12);
    val_used  = pu .* Ku_e;
    ushare    = val_used ./ max(Kn_e + val_used, 1e-12);

    wq = @(x, w) sum(x(:).*w(:)) / max(sum(w(:)), 1e-300);

    CUTFILE = fullfile(fileparts(sd), 'common_cuts.mat');
    if e == 1
        mtot_c  = me_raw + mw_raw;
        [~, oc] = sort(a_s);
        cmc     = cumsum(mtot_c(oc)) / max(sum(mtot_c), 1e-300);
        edges_a = zeros(1,3);
        for q = 1:3
            idx = find(cmc >= q/4, 1, 'first');
            if isempty(idx), idx = numel(oc); end
            edges_a(q) = a_s(oc(idx));
        end
        imed_c = find(cmc >= 0.5, 1, 'first');
        if isempty(imed_c), imed_c = numel(oc); end
        a_med_common = a_s(oc(imed_c));
        [~, oe] = sort(a_s);
        cme_c   = cumsum(me(oe)) / max(sum(me), 1e-300);
        edges_d_common = zeros(1,9);
        for q = 1:9
            idx = find(cme_c >= q/10, 1, 'first');
            if isempty(idx), idx = numel(oe); end
            edges_d_common(q) = a_s(oe(idx));
        end
        save(CUTFILE, 'edges_a', 'a_med_common', 'edges_d_common');
        fprintf('\n  [cuts] COMMON cut points written to %s\n', CUTFILE);
        fprintf('         wealth edges = [%.6g  %.6g  %.6g]\n', edges_a);
        fprintf('         median wealth = %.6g   (the "poor" threshold)\n', a_med_common);
    else
        Cc = load(CUTFILE);
        edges_a = Cc.edges_a;  a_med_common = Cc.a_med_common;
        edges_d_common = Cc.edges_d_common;
    end
    aq = ones(size(a_s));
    aq(a_s > edges_a(1)) = 2;
    aq(a_s > edges_a(2)) = 3;
    aq(a_s > edges_a(3)) = 4;

    zt = ones(size(z_s));
    zt(z_s > round(k/3))   = 2;
    zt(z_s > round(2*k/3)) = 3;

    tabA    = zeros(4,3);
    tabENT  = zeros(4,3);
    tabPOP  = zeros(4,3);
    tabMASSE= zeros(4,3);
    tabZ    = zeros(4,3);
    tabMRPK = zeros(4,3);
    tabMU   = zeros(4,3);
    for q = 1:4
        for t = 1:3
            sel = (aq == q) & (zt == t);
            tabA(q,t)     = wq(ushare(sel), me(sel));
            me_c          = sum(me_pop(sel));
            mw_c          = sum(mw_pop(sel));
            tabMASSE(q,t) = me_c;
            tabPOP(q,t)   = me_c + mw_c;
            tabENT(q,t)   = me_c / max(me_c + mw_c, 1e-300);
            tabZ(q,t)     = wq(Eeff(sel),  me(sel));
            tabMRPK(q,t)  = wq(MRPK(sel),  me(sel));
            tabMU(q,t)    = wq(mu_e(sel),  me(sel));
        end
    end

    edges_d = edges_d_common;
    adec = ones(size(a_s));
    for q = 1:9
        adec(a_s > edges_d(q)) = q + 1;
    end
    tabMU10 = nan(10,3);
    tabU10  = nan(10,3);
    tabM10  = zeros(10,3);
    for dd = 1:10
        for t = 1:3
            sel = (adec == dd) & (zt == t);
            tabM10(dd,t) = sum(me(sel));
            if tabM10(dd,t) > 1e-14
                tabMU10(dd,t) = wq(mu_e(sel),   me(sel));
                tabU10(dd,t)  = wq(ushare(sel), me(sel));
            end
        end
    end

    qcut = @(x) local_quartile(x, me);
    Lq = qcut(L_e);
    bindQ = zeros(4,1);
    empQ  = zeros(4,1);
    for q = 1:4
        sel = (Lq == q);
        bindQ(q) = wq(double(bind_e(sel)), me(sel));
        empQ(q)  = wq(L_e(sel), me(sel));
    end

    Zq = qcut(Eeff);
    Kq = zeros(4,1);
    Yq = zeros(4,1);
    for q = 1:4
        Kq(q) = sum(me_raw(Zq == q) .* Kagg_e(Zq == q));
        Yq(q) = sum(me_raw(Zq == q) .* Y_e(Zq == q));
    end

    mtot = me_raw + mw_raw;
    a_med = a_med_common;
    cz    = cumsum(Perg(:)) / max(sum(Perg), 1e-300);
    topz  = find(cz > 0.90);
    selPT = (a_s <= a_med) & ismember(z_s, topz);
    selT  = ismember(z_s, topz);
    pt_rate  = sum(me_raw(selPT)) / max(sum(mtot(selPT)), 1e-300);
    t10_rate = sum(me_raw(selT))  / max(sum(mtot(selT)),  1e-300);
    rich_pt  = (a_s > a_med) & ismember(z_s, topz);
    ptr_rich = sum(me_raw(rich_pt)) / max(sum(mtot(rich_pt)), 1e-300);

    Ktot = sum(me_raw .* Kagg_e);
    Kn_tot = sum(me_raw .* Kn_e);
    Ku_tot = sum(me_raw .* Ku_e);
    if Kn_tot > 0 && Ku_tot > 0
        rho_pool = Ku_tot / Kn_tot;
        Psi_pool = ( phi^(1/gamma) + (1-phi)^(1/gamma) ...
            * rho_pool^((gamma-1)/gamma) )^(gamma/(gamma-1));
        Ktot_pool = Psi_pool * Kn_tot;
    else
        Ktot_pool = Ktot;
    end
    ces_slack = 100 * log(max(Ktot_pool,1e-300) / max(Ktot,1e-300));
    Ltot = sum(me_raw .* L_e);
    Yact = sum(me_raw .* Y_e);
    Sz   = sum(me_raw .* Eeff.^(1/(1-eta)));
    Yeff = (Ltot^alpha * Ktot^(1-alpha))^eta * Sz^(1-eta);
    intensive_gain = 100 * log(max(Yeff,1e-300) / max(Yact,1e-300));
    Yeff_pool = (Ltot^alpha * Ktot_pool^(1-alpha))^eta * Sz^(1-eta);
    intensive_gain_pool = 100 * log(max(Yeff_pool,1e-300) / max(Yact,1e-300));

    AK      = size(INV_AGE_COH, 2);
    age_u   = nan(AK,1);
    age_mu  = nan(AK,1);
    age_m   = nan(AK,1);
    age_sdm = nan(AK,1);
    age_mrpk= nan(AK,1);
    lm      = log(max(MRPK, 1e-12));
    for aa = 1:AK
        c = INV_AGE_COH(:,aa);
        age_m(aa) = sum(c);
        if age_m(aa) > 1e-14
            age_u(aa)   = wq(ushare, c);
            age_mu(aa)  = wq(mu_e,   c);
            age_mrpk(aa)= wq(MRPK,   c);
            wgt  = c / sum(c);
            mbar = sum(wgt .* lm);
            age_sdm(aa) = sqrt(max(sum(wgt .* (lm - mbar).^2), 0));
        end
    end
    age_m = age_m / max(sum(INV_AGE_MASS), 1e-300);

    mbar_all = wq(lm, me);
    sd_all   = sqrt(max(wq((lm - mbar_all).^2, me), 0));

    R(e).name   = ECON{e,1};
    R(e).ushare = ushare;  R(e).mu = mu_e;  R(e).MRPK = MRPK;
    R(e).me = me;  R(e).a = a_s;  R(e).z = z_s;  R(e).aq = aq;  R(e).zt = zt;
    R(e).age_used = age_u;  R(e).age_mu = age_mu;  R(e).age_mass = age_m;
    R(e).age_mrpk = age_mrpk;  R(e).age_sdm = age_sdm;  R(e).tabA = tabA;
    R(e).tabENT = tabENT;  R(e).tabPOP = tabPOP;  R(e).tabMASSE = tabMASSE;
    R(e).tabZ = tabZ;  R(e).tabMRPK = tabMRPK;  R(e).tabMU = tabMU;
    R(e).tabMU10 = tabMU10;  R(e).tabU10 = tabU10;  R(e).tabM10 = tabM10;
    R(e).Ktot = Ktot;  R(e).Ltot = Ltot;  R(e).Yact = Yact;  R(e).Yeff = Yeff;
    R(e).intensive_gain = intensive_gain;
    R(e).Sz = Sz;  R(e).TFPeff_set = Sz^(1-eta);
    R(e).bindQ = bindQ;  R(e).empQ = empQ;  R(e).Kq = Kq;  R(e).Yq = Yq;
    R(e).pt_rate = pt_rate;  R(e).t10_rate = t10_rate;
    R(e).ptr_rich = ptr_rich;  R(e).a_med = a_med;
    R(e).sd_logMRPK = sd_all;
    R(e).W = W;  R(e).r = r;  R(e).Rnew = Rnew;  R(e).qu = qu;  R(e).pu = pu;
    R(e).theta = theta;

    AG = aggregates();
    R(e).Y      = AG.Y;
    R(e).TFP    = AG.TFP;
    R(e).Kagg   = AG.Kagg;
    R(e).constr = AG.constr;
    R(e).ent    = AG.ent;

    gb = gaps_calibration_last;
    R(e).gap_Knew  = 100*gb(1);
    R(e).gap_Labor = 100*gb(2);
    fprintf('  [%s] Y=%.5f  TFP=%.5f  Kagg=%.4f  ent=%.4f  constrained=%.1f%%\n', ...
        ECON{e,1}, AG.Y, AG.TFP, AG.Kagg, AG.ent, 100*AG.constr);
    fprintf('  [%s] market gaps at these prices: Knew %+.3f%%  Labor %+.3f%%\n', ...
        ECON{e,1}, R(e).gap_Knew, R(e).gap_Labor);
    if max(abs([R(e).gap_Knew, R(e).gap_Labor])) > 1.0
        warning(['analysis_used_capital: the %s economy is NOT cleared at these prices ' ...
            '(gaps Knew %+.2f%%, Labor %+.2f%%). The statistics below are off-equilibrium.'], ...
            ECON{e,1}, R(e).gap_Knew, R(e).gap_Labor);
    end
end

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TABLE A -- used-capital VALUE share, wealth quartile x ability tercile\n');
fprintf('  entry = pu*Ku / (Kn + pu*Ku), entrepreneur-mass weighted\n');
fprintf('================================================================\n');
for e = 1:NE
    fprintf('\n  %s\n', upper(R(e).name));
    fprintf('  %-18s %12s %12s %12s\n', 'wealth quartile', 'low abil', 'mid abil', 'high abil');
    for q = 1:4
        fprintf('  Q%-17d', q);
        for t = 1:3
            fprintf('%12.4f', R(e).tabA(q,t));
        end
        fprintf('\n');
    end
end
fprintf('\n  Q1 = poorest.  A high entry means that cell finances itself from the\n');
fprintf('  second-hand market.  If the poor-and-talented cell (Q1, high abil) is the\n');
fprintf('  largest, used capital is doing the entry work the thesis claims.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TABLE A2 -- ENTRY RATE by wealth quartile x ability tercile\n');
fprintf('  entry rate = entrepreneurs / (entrepreneurs + workers) within the cell\n');
fprintf('================================================================\n');
for e = 1:NE
    fprintf('\n  %s -- entry rate\n', upper(R(e).name));
    fprintf('  %-18s %12s %12s %12s %12s\n', 'wealth quartile', 'low abil', 'mid abil', 'high abil', 'cell pop');
    for q = 1:4
        fprintf('  Q%-17d', q);
        for t = 1:3
            fprintf('%12.4f', R(e).tabENT(q,t));
        end
        fprintf('%12.4f\n', sum(R(e).tabPOP(q,:)));
    end
end
fprintf('\n  This is the selection margin in full. The benchmark-minus-CF difference\n');
fprintf('  in this matrix identifies exactly WHICH cells used capital lets in.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TFP DECOMPOSITION -- intensive (allocation) vs extensive (selection)\n');
fprintf('  y = z*(L^alpha K^(1-alpha))^eta, so with the operator set, K and L HELD\n');
fprintf('  FIXED the efficient allocation gives every firm scale proportional to\n');
fprintf('  z^(1/(1-eta)) and Y_eff = (L^alpha K^(1-alpha))^eta * [sum m z^(1/(1-eta))]^(1-eta).\n');
fprintf('================================================================\n\n');
fprintf('  %-46s %14s %14s\n', 'quantity', 'friction', 'frictionless');
fprintf('  %-46s %14.5f %14.5f\n', 'Y actual', R(1).Yact, R(2).Yact);
fprintf('  %-46s %14.5f %14.5f\n', 'Y if capital reshuffled among SAME operators', R(1).Yeff, R(2).Yeff);
fprintf('  %-46s %13.2f%% %13.2f%%\n', 'INTENSIVE gain available (log pts)', ...
    R(1).intensive_gain, R(2).intensive_gain);
tot_loss = 100 * log(R(2).TFP / R(1).TFP);
fprintf('\n  %-46s %14.2f\n', 'TOTAL loss, friction -> frictionless (log pts)', tot_loss);
fprintf('  %-46s %14.2f\n', '  of which INTENSIVE (reallocation)', R(1).intensive_gain);
fprintf('  %-46s %14.2f\n', '  of which EXTENSIVE (selection)     ', tot_loss - R(1).intensive_gain);
fprintf('\n  The intensive term is what a planner could gain by reshuffling capital among\n');
fprintf('  the firms that already operate. The remainder is what only a change in WHO\n');
fprintf('  operates, and the price movements that follow, can deliver.\n');
fprintf('  The frictionless column should show an intensive gain of ~0: with no\n');
fprintf('  constraint, capital is already efficiently allocated. Treat a large value\n');
fprintf('  there as a diagnostic failure, not a result.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TABLE B -- by firm age: used share, shadow cost mu, MRPK, dispersion\n');
fprintf('================================================================\n');
for e = 1:NE
    fprintf('\n  %s\n', upper(R(e).name));
    fprintf('  %5s %11s %11s %10s %12s %12s\n', 'age', 'mass', 'used share', 'mean mu', 'mean MRPK', 'sd log MRPK');
    for aa = 1:min(size(INV_AGE_COH,2), 20)
        if isfinite(R(e).age_used(aa))
            fprintf('  %5d %11.5f %11.4f %10.4f %12.4f %12.4f\n', ...
                aa, R(e).age_mass(aa), R(e).age_used(aa), R(e).age_mu(aa), ...
                R(e).age_mrpk(aa), R(e).age_sdm(aa));
        end
    end
end

AK   = size(INV_AGE_COH, 2);
yidx = 1:min(5,AK);
oidx = min(15,AK):min(AK,40);

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SELECTION AND ALLOCATION -- friction against frictionless\n');
fprintf('================================================================\n\n');
fprintf('  %-42s %14s %14s\n', 'statistic', 'friction', 'frictionless');
fprintf('  %-42s %14.5f %14.5f\n', 'output Y', R(1).Y, R(2).Y);
fprintf('  %-42s %14.5f %14.5f\n', 'TFP', R(1).TFP, R(2).TFP);
fprintf('  %-42s %14.4f %14.4f\n', 'entrepreneur share', R(1).ent, R(2).ent);
fprintf('  %-42s %13.1f%% %13.1f%%\n', 'constrained entrepreneurs', 100*R(1).constr, 100*R(2).constr);
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), all firms', R(1).sd_logMRPK, R(2).sd_logMRPK);
for e = 1:NE
    my(e) = mean(R(e).age_mrpk(yidx),'omitnan');
    mo(e) = mean(R(e).age_mrpk(oidx),'omitnan');
    sy(e) = mean(R(e).age_sdm(yidx),'omitnan');
    so(e) = mean(R(e).age_sdm(oidx),'omitnan');
end
fprintf('  %-42s %14.4f %14.4f\n', 'mean MRPK, ages 1-5',  my(1), my(2));
fprintf('  %-42s %14.4f %14.4f\n', 'mean MRPK, ages 15+',  mo(1), mo(2));
fprintf('  %-42s %+13.1f%% %+13.1f%%\n', 'young MRPK premium (Midrigan-Xu: +73%)', ...
    100*(my(1)/mo(1)-1), 100*(my(2)/mo(2)-1));
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), ages 1-5', sy(1), sy(2));
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), ages 15+', so(1), so(2));
fprintf('  %-42s %14.4f %14.4f\n', 'mean mu, ages 1-5', mean(R(1).age_mu(yidx),'omitnan'), mean(R(2).age_mu(yidx),'omitnan'));
fprintf('  %-42s %14.4f %14.4f\n', 'used share, ages 1-5', mean(R(1).age_used(yidx),'omitnan'), mean(R(2).age_used(yidx),'omitnan'));
fprintf('\n  The frictionless column is the benchmark of perfect allocation: mu = 0 and\n');
fprintf('  sd(log MRPK) collapses to whatever span-of-control curvature alone implies.\n');
fprintf('  The gap between the columns IS the friction, measured age by age.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  CHEN (2023) DIAGNOSTICS -- computed WITHIN this economy\n');
fprintf('================================================================\n');

mrpk_star = R(2).age_mrpk(min(15,AK));
if ~isfinite(mrpk_star), mrpk_star = mean(R(2).age_mrpk,'omitnan'); end
tolsf = 0.05;
psf = find(R(1).age_mrpk <= mrpk_star*(1+tolsf) & isfinite(R(1).age_mrpk), 1, 'first');
fprintf('\n  PERIODS TO SELF-FINANCE\n');
fprintf('    benchmark of an unconstrained firm: MRPK = %.4f (frictionless twin)\n', mrpk_star);
fprintf('    first age at which mean MRPK falls within %.0f%% of it: ', 100*tolsf);
if isempty(psf)
    fprintf('NOT within %d years\n', AK);
    fprintf('    gap still %+.1f%% at age %d -- these firms never fully self-finance.\n', ...
        100*(R(1).age_mrpk(min(20,AK))/mrpk_star - 1), min(20,AK));
else
    fprintf('%d\n', psf);
end
fprintf('    (Chen reports 35 with managerial inputs, 19 without.)\n');

fprintf('\n  PCT FINANCIALLY CONSTRAINED, by employment quartile\n');
fprintf('    %-10s %12s %12s %12s\n', '', 'friction', 'frictionless', 'mean L');
for q = 1:4
    fprintf('    Q%-9d %12.3f %12.3f %12.3f\n', q, R(1).bindQ(q), R(2).bindQ(q), R(1).empQ(q));
end
fprintf('    Q1 = smallest by employment.  Chen: 0.309, 0.423, 0.473, 0.593.\n');

fprintf('\n  CAPITAL INPUT RELATIVE TO BASELINE, by firm-TFP quartile\n');
fprintf('    (what removing the constraint does to each quartile''s capital)\n');
fprintf('    %-10s %14s %14s %14s\n', '', 'K friction', 'K frictionless', 'relative');
for q = 1:4
    fprintf('    Q%-9d %14.5f %14.5f %+13.3f\n', q, R(1).Kq(q), R(2).Kq(q), ...
        R(2).Kq(q)/max(R(1).Kq(q),1e-12) - 1);
end
fprintf('    Q1 = lowest TFP.  Chen: -0.300, -0.269, -0.134, +0.216: capital leaves the\n');
fprintf('    unproductive and goes to the productive. That IS the reallocation.\n');

fprintf('\n  SELECTION INDEX of the operating set (TFP under perfect allocation)\n');
fprintf('    friction     sum m z^(1/(1-eta)) = %.5f   -> TFP_eff = %.5f\n', R(1).Sz, R(1).TFPeff_set);
fprintf('    frictionless sum m z^(1/(1-eta)) = %.5f   -> TFP_eff = %.5f\n', R(2).Sz, R(2).TFPeff_set);

fprintf('\n  ENTRY OF THE POOR AND TALENTED -- fixed cut, comparable across economies\n');
fprintf('  poor = below the POPULATION median wealth; talented = top 10%% of the\n');
fprintf('  ergodic ability distribution. Neither cut moves between economies.\n\n');
fprintf('    %-42s %13s %13s\n', 'statistic', 'friction', 'frictionless');
fprintf('    %-42s %13.4f %13.4f\n', 'poor & talented entry rate',  R(1).pt_rate,  R(2).pt_rate);
fprintf('    %-42s %13.4f %13.4f\n', 'RICH & talented entry rate',  R(1).ptr_rich, R(2).ptr_rich);
fprintf('    %-42s %13.4f %13.4f\n', 'all talented entry rate',     R(1).t10_rate, R(2).t10_rate);
fprintf('    %-42s %13.4f %13.4f\n', 'wealth gap (rich - poor)', ...
    R(1).ptr_rich-R(1).pt_rate, R(2).ptr_rich-R(2).pt_rate);
fprintf('\n    The frictionless column is the no-distortion benchmark: with unlimited\n');
fprintf('    borrowing, wealth should barely matter for a talented agent, so the\n');
fprintf('    wealth gap should collapse. The friction column shows how far the\n');
fprintf('    collateral constraint pushes entry away from that.\n');
fprintf('    Poor talented agents enter at %.1f%% of the rate perfect credit would give.\n', ...
    100*R(1).pt_rate/max(R(2).pt_rate,1e-12));
fprintf('    This is the TFP the operating set WOULD deliver if capital were allocated\n');
fprintf('    efficiently, so it isolates who operates from how capital is spread.\n');

fprintf('\n');
fprintf('================================================================\n');
fprintf('  MISALLOCATION LOSS  100*log(TFP frictionless / TFP friction)\n');
fprintf('  Midrigan convention. A WITHIN-economy ratio, so it is free of level effects.\n');
fprintf('  Both economies here carry used capital; the only difference is the\n');
fprintf('  collateral constraint. So this is the loss the CONSTRAINT causes, given\n');
fprintf('  that a second-hand market already exists.\n');
fprintf('================================================================\n\n');

loss_saved = 100 * log(NFR.frictionless.TFP / NFR.base.TFP);
loss_self  = 100 * log(R(2).TFP / R(1).TFP);

fprintf('  AUTHORITATIVE -- from nofriction.m, where BOTH economies were GE-cleared:\n');
fprintf('    TFP with the constraint     = %.6f   (theta = %.4f, W=%.4f r=%.4f)\n', ...
    NFR.base.TFP, NFR.base.theta, NFR.base.W, NFR.base.r);
fprintf('    TFP without the constraint  = %.6f   (theta = %.4g, W=%.4f r=%.4f)\n', ...
    NFR.frictionless.TFP, NFR.frictionless.theta, NFR.frictionless.W, NFR.frictionless.r);
fprintf('    misallocation loss          = %.4f log points\n', loss_saved);
fprintf('    output gain from perfect credit = %+.2f%%\n', NFR.gain_Y_pct);
fprintf('    residual clearing gaps: base Knew %+.2f%% Labor %+.2f%%\n', ...
    NFR.base.gap_Knew, NFR.base.gap_Labor);

fprintf('\n  CHECK -- recomputed here at those same cleared prices:\n');
fprintf('    TFP %.6f and %.6f  ->  loss %.4f   (difference %.2e)\n', ...
    R(1).TFP, R(2).TFP, loss_self, loss_self - loss_saved);
if abs(loss_self - loss_saved) > 0.05
    fprintf('  ^ THESE DISAGREE. This analysis and the GE-cleared run are not the same\n');
    fprintf('    economy. Trust the authoritative block and re-run nofriction.m.\n');
else
    fprintf('    Agreement confirms the analysis sits on the GE-cleared equilibrium.\n');
end

A = struct('econ',{R}, 'misalloc_loss',loss_self, 'misalloc_saved',loss_saved, 'theta_cal',theta_cal);
save('analysis_used_capital_result.mat','A');

f = figure('Visible','off','Position',[100 100 1200 380]);
nplot = min(size(INV_AGE_COH,2), 25);

subplot(1,3,1);
yyaxis left
plot(1:nplot, 100*R(1).age_used(1:nplot), '-o', 'LineWidth', 1.6, 'MarkerSize', 4);
ylabel('used-capital value share (%)');
yyaxis right
plot(1:nplot, R(1).age_mu(1:nplot), '-s', 'LineWidth', 1.6, 'MarkerSize', 4);
ylabel('shadow cost \mu');
xlabel('firm age (years since entry)');
title('Friction: used capital and \mu over firm age');
grid on;

subplot(1,3,2);
plot(1:nplot, R(1).age_sdm(1:nplot), '-o', 'LineWidth', 1.6, 'MarkerSize', 4); hold on;
plot(1:nplot, R(2).age_sdm(1:nplot), '-s', 'LineWidth', 1.6, 'MarkerSize', 4);
xlabel('firm age (years since entry)');
ylabel('sd(log MRPK)');
legend({'friction','frictionless'}, 'Location','northeast');
title('Allocation: MRPK dispersion by age');
grid on;

subplot(1,3,3);
b = bar(100*R(1).tabA);
set(gca,'XTickLabel',{'Q1 poorest','Q2','Q3','Q4 richest'});
ylabel('used-capital value share (%)');
legend({'low ability','mid ability','high ability'},'Location','northeast');
title('Who uses used capital: wealth x ability');
grid on;

print(f, 'fig_used_capital.png', '-dpng', '-r150');
close(f);

fprintf('\n  [saved] analysis_used_capital_result.mat\n');
fprintf('  [saved] fig_used_capital.png\n');

function q = local_quartile(x, w)
[~, ord] = sort(x(:));
cm  = cumsum(w(ord)) / max(sum(w), 1e-300);
q   = ones(size(x(:)));
for j = 1:3
    idx = find(cm >= j/4, 1, 'first');
    if isempty(idx), idx = numel(ord); end
    q(x(:) > x(ord(idx))) = j + 1;
end
q = reshape(q, size(x));
end
