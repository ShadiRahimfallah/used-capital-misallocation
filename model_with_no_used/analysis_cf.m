clear; clear global; clc;

global rho alpha eta beta delta W r smin smax fspacee fspacew fspace ce cw ...
    cxe cxw cse csw egrid P smine sminw smaxw theta xai F Lbar phi gamma se ...
    Rnew eps_model productivity_weight Perg k exit_prob smaxe qu zeta kappa ...
    pu use_calib_sim_size calib_mode DtoY ygrid Py Ny rho_y sigma_y Ey piy ...
    sigma_crra lambda_u sigma_perm eta_p omega fast_diagnostic ...
    wage_gini_target firm_exit_target top10_emp_share_target ...
    xsim_prev_for_amax amax_calib_mult hybrid_use_nested_clearing ...
    hybrid_inner_clearing_skip hybrid_nested_price_lb hybrid_nested_price_ub ...
    gaps_calibration_last LOCKOUT_EXOG_DEATH nested_clearing_quiet ...
    nested_clearing_max_iter nested_clearing_tol nested_clearing_tol_Knew ...
    calib_moment_weights USE_INVARIANT_SEED SIM_NF SIM_TINIT LAST_MOMENTS ...
    psi_u WGINI_USE_TRIM INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT ...
    INV_ME INV_MW INV_AG COMPUTE_AGE INV_AGE_MASS INV_AGE_COH INV_AGE_KEEP

sd = fileparts(mfilename('fullpath'));
cd(sd);
startup;

parameters_cf_calib;

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
WGINI_USE_TRIM             = false;
hybrid_nested_price_lb     = [0.35; 0.005; 0.04];
hybrid_nested_price_ub     = [2.00; 0.99*(1/beta-1); 0.26];

COMPUTE_AGE   = true;
INV_AGE_KEEP  = 40;

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic);
[ygrid, Py, Ey, piy]          = rouwenhorst_y(Ny, rho_y, sigma_y);

theta_cal = theta;
W_cal     = W;
r_cal     = r;
Rnew_cal  = Rnew;

NFFILE = fullfile(sd, 'nofriction_result.mat');
if ~exist(NFFILE, 'file')
    error(['analysis_cf: %s not found. Run nofriction.m in this folder first -- this ' ...
        'analysis compares the CF against its own GE-cleared frictionless twin.'], NFFILE);
end
NF  = load(NFFILE);
fn  = fieldnames(NF);
NFR = NF.(fn{1});
if abs(NFR.base.theta - theta_cal) > 1e-6
    warning(['analysis_cf: nofriction_result.mat was built at theta=%.6f but ' ...
        'parameters_cf_calib has theta=%.6f. Not the same economy.'], NFR.base.theta, theta_cal);
end

ECON = { 'friction',     theta_cal, W_cal,              r_cal,              Rnew_cal
         'frictionless', 1e5,       NFR.frictionless.W, NFR.frictionless.r, NFR.frictionless.Rnew };
NE = size(ECON, 1);

fprintf('================================================================\n');
fprintf('  NO-USED-CAPITAL CF vs ITS FRICTIONLESS TWIN -- selection and allocation\n');
fprintf('  theta_cal=%.4f eta_p=%.4f omega=%.4f sigma_y=%.4f beta=%.4f\n', ...
    theta_cal, eta_p, omega, sigma_y, beta);
fprintf('  Used capital is ABSENT in both: kappa=0, lambda_u=0, phi=1, qu=pu=0.\n');
for i = 1:NE
    fprintf('  %-13s theta=%10.4g  W=%.4f r=%.4f Rnew=%.4f\n', ...
        ECON{i,1}, ECON{i,2}, ECON{i,3}, ECON{i,4}, ECON{i,5});
end
fprintf('  Prices held at each economy''s own GE-cleared equilibrium; nothing is re-cleared.\n');
fprintf('================================================================\n');

R = struct('name',{},'mu',{},'MRPK',{},'me',{},'a',{},'z',{},'aq',{},'zt',{}, ...
    'tabMU',{},'tabBIND',{},'age_mu',{},'age_mass',{},'age_mrpk',{},'age_sdm',{}, ...
    'sd_logMRPK',{},'ent',{},'Y',{},'TFP',{},'Kagg',{},'constr',{}, ...
    'W',{},'r',{},'Rnew',{},'theta',{},'gap_Knew',{},'gap_Labor',{});

for e = 1:NE
    theta = ECON{e,2};
    W     = ECON{e,3};
    r     = ECON{e,4};
    Rnew  = ECON{e,5};
    qu    = 0;
    pu    = 0;

    fprintf('\n\n----------------------------------------------------------------\n');
    fprintf('  SOLVING: %s   (theta=%.4g, W=%.4f, r=%.4f)\n', ECON{e,1}, theta, W, r);
    fprintf('----------------------------------------------------------------\n');

    x_full = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    for w = 1:3
        evaluate_economy(x_full);
    end
    evaluate_economy(x_full);

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

    MRPK = (1 - alpha) * eta * Y_e ./ max(Kagg_e, 1e-12);

    wq = @(x, w) sum(x(:).*w(:)) / max(sum(w(:)), 1e-300);

    CUTFILE = fullfile(fileparts(sd), 'common_cuts.mat');
    if ~exist(CUTFILE, 'file')
        error(['analysis_cf (no used capital): %s not found. Run the BENCHMARK analysis ' ...
            'FIRST -- it writes the common cut points.'], CUTFILE);
    end
    Cc = load(CUTFILE);
    edges_a = Cc.edges_a;  a_med_common = Cc.a_med_common;
    edges_d_common = Cc.edges_d_common;
    aq = ones(size(a_s));
    aq(a_s > edges_a(1)) = 2;
    aq(a_s > edges_a(2)) = 3;
    aq(a_s > edges_a(3)) = 4;

    zt = ones(size(z_s));
    zt(z_s > round(k/3))   = 2;
    zt(z_s > round(2*k/3)) = 3;

    tabMU    = zeros(4,3);
    tabBIND  = zeros(4,3);
    tabENT   = zeros(4,3);
    tabPOP   = zeros(4,3);
    tabMASSE = zeros(4,3);
    tabZ     = zeros(4,3);
    tabMRPK  = zeros(4,3);
    for q = 1:4
        for t = 1:3
            sel = (aq == q) & (zt == t);
            tabMU(q,t)   = wq(mu_e(sel),   me(sel));
            tabBIND(q,t) = wq(double(bind_e(sel)), me(sel));
            me_c          = sum(me_pop(sel));
            mw_c          = sum(mw_pop(sel));
            tabMASSE(q,t) = me_c;
            tabPOP(q,t)   = me_c + mw_c;
            tabENT(q,t)   = me_c / max(me_c + mw_c, 1e-300);
            tabZ(q,t)     = wq(Eeff(sel), me(sel));
            tabMRPK(q,t)  = wq(MRPK(sel), me(sel));
        end
    end

    edges_d = edges_d_common;
    adec = ones(size(a_s));
    for q = 1:9
        adec(a_s > edges_d(q)) = q + 1;
    end
    tabMU10 = nan(10,3);
    tabM10  = zeros(10,3);
    for dd = 1:10
        for t = 1:3
            sel = (adec == dd) & (zt == t);
            tabM10(dd,t) = sum(me(sel));
            if tabM10(dd,t) > 1e-14
                tabMU10(dd,t) = wq(mu_e(sel), me(sel));
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
    Ltot = sum(me_raw .* L_e);
    Yact = sum(me_raw .* Y_e);
    Sz   = sum(me_raw .* Eeff.^(1/(1-eta)));
    Yeff = (Ltot^alpha * Ktot^(1-alpha))^eta * Sz^(1-eta);
    intensive_gain = 100 * log(max(Yeff,1e-300) / max(Yact,1e-300));

    AK       = size(INV_AGE_COH, 2);
    age_mu   = nan(AK,1);
    age_m    = nan(AK,1);
    age_sdm  = nan(AK,1);
    age_mrpk = nan(AK,1);
    lm       = log(max(MRPK, 1e-12));
    for aa = 1:AK
        c = INV_AGE_COH(:,aa);
        age_m(aa) = sum(c);
        if age_m(aa) > 1e-14
            age_mu(aa)   = wq(mu_e, c);
            age_mrpk(aa) = wq(MRPK, c);
            wgt  = c / sum(c);
            mbar = sum(wgt .* lm);
            age_sdm(aa) = sqrt(max(sum(wgt .* (lm - mbar).^2), 0));
        end
    end
    age_m = age_m / max(sum(INV_AGE_MASS), 1e-300);

    mbar_all = wq(lm, me);
    sd_all   = sqrt(max(wq((lm - mbar_all).^2, me), 0));

    AG = aggregates();

    R(e).name = ECON{e,1};
    R(e).mu = mu_e;  R(e).MRPK = MRPK;  R(e).me = me;
    R(e).a = a_s;  R(e).z = z_s;  R(e).aq = aq;  R(e).zt = zt;
    R(e).tabMU = tabMU;  R(e).tabBIND = tabBIND;
    R(e).tabENT = tabENT;  R(e).tabPOP = tabPOP;  R(e).tabMASSE = tabMASSE;
    R(e).tabZ = tabZ;  R(e).tabMRPK = tabMRPK;
    R(e).tabMU10 = tabMU10;  R(e).tabM10 = tabM10;
    R(e).Ktot = Ktot;  R(e).Ltot = Ltot;  R(e).Yact = Yact;  R(e).Yeff = Yeff;
    R(e).intensive_gain = intensive_gain;
    R(e).Sz = Sz;  R(e).TFPeff_set = Sz^(1-eta);
    R(e).bindQ = bindQ;  R(e).empQ = empQ;  R(e).Kq = Kq;  R(e).Yq = Yq;
    R(e).pt_rate = pt_rate;  R(e).t10_rate = t10_rate;
    R(e).ptr_rich = ptr_rich;  R(e).a_med = a_med;
    R(e).age_mu = age_mu;  R(e).age_mass = age_m;
    R(e).age_mrpk = age_mrpk;  R(e).age_sdm = age_sdm;
    R(e).sd_logMRPK = sd_all;
    R(e).Y = AG.Y;  R(e).TFP = AG.TFP;  R(e).Kagg = AG.Kagg;
    R(e).constr = AG.constr;  R(e).ent = AG.ent;
    R(e).W = W;  R(e).r = r;  R(e).Rnew = Rnew;  R(e).theta = theta;

    gb = gaps_calibration_last;
    R(e).gap_Knew  = 100*gb(1);
    R(e).gap_Labor = 100*gb(2);
    fprintf('  [%s] Y=%.5f  TFP=%.5f  Kagg=%.4f  ent=%.4f  constrained=%.1f%%\n', ...
        ECON{e,1}, AG.Y, AG.TFP, AG.Kagg, AG.ent, 100*AG.constr);
    fprintf('  [%s] market gaps at these prices: Knew %+.3f%%  Labor %+.3f%%\n', ...
        ECON{e,1}, R(e).gap_Knew, R(e).gap_Labor);
    if max(abs([R(e).gap_Knew, R(e).gap_Labor])) > 1.0
        warning(['analysis_cf: the %s economy is NOT cleared at these prices ' ...
            '(Knew %+.2f%%, Labor %+.2f%%).'], ECON{e,1}, R(e).gap_Knew, R(e).gap_Labor);
    end
end

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TABLE A -- the constraint by wealth quartile x ability tercile\n');
fprintf('  There is no used-capital share to report here: the CF has no second-hand\n');
fprintf('  market, so the comparable object is how hard the constraint bites.\n');
fprintf('================================================================\n');
for e = 1:NE
    fprintf('\n  %s -- mean shadow cost mu\n', upper(R(e).name));
    fprintf('  %-18s %12s %12s %12s\n', 'wealth quartile', 'low abil', 'mid abil', 'high abil');
    for q = 1:4
        fprintf('  Q%-17d', q);
        for t = 1:3
            fprintf('%12.4f', R(e).tabMU(q,t));
        end
        fprintf('\n');
    end
    fprintf('\n  %s -- share with a BINDING constraint\n', upper(R(e).name));
    fprintf('  %-18s %12s %12s %12s\n', 'wealth quartile', 'low abil', 'mid abil', 'high abil');
    for q = 1:4
        fprintf('  Q%-17d', q);
        for t = 1:3
            fprintf('%12.4f', R(e).tabBIND(q,t));
        end
        fprintf('\n');
    end
end

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
fprintf('\n  Compare cell by cell against the benchmark''s Table A2 to see exactly which\n');
fprintf('  agents used capital lets into entrepreneurship.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TFP DECOMPOSITION -- intensive (allocation) vs extensive (selection)\n');
fprintf('  Y_eff = (L^alpha K^(1-alpha))^eta * [sum m z^(1/(1-eta))]^(1-eta), the output\n');
fprintf('  attainable by reshuffling capital among the SAME operators.\n');
fprintf('================================================================\n\n');
fprintf('  %-46s %14s %14s\n', 'quantity', 'friction', 'frictionless');
fprintf('  %-46s %14.5f %14.5f\n', 'Y actual', R(1).Yact, R(2).Yact);
fprintf('  %-46s %14.5f %14.5f\n', 'Y if capital reshuffled among SAME operators', R(1).Yeff, R(2).Yeff);
fprintf('  %-46s %13.2f%% %13.2f%%\n', 'INTENSIVE gain available (log pts)', ...
    R(1).intensive_gain, R(2).intensive_gain);
tot_loss_cf = 100 * log(R(2).TFP / R(1).TFP);
fprintf('\n  %-46s %14.2f\n', 'TOTAL loss, friction -> frictionless (log pts)', tot_loss_cf);
fprintf('  %-46s %14.2f\n', '  of which INTENSIVE (reallocation)', R(1).intensive_gain);
fprintf('  %-46s %14.2f\n', '  of which EXTENSIVE (selection) + GE', tot_loss_cf - R(1).intensive_gain);
fprintf('\n  The frictionless column should show an intensive gain of ~0. A large value\n');
fprintf('  there is a diagnostic failure, not a result.\n');

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  TABLE B -- by firm age: shadow cost, MRPK, dispersion\n');
fprintf('================================================================\n');
for e = 1:NE
    fprintf('\n  %s\n', upper(R(e).name));
    fprintf('  %5s %11s %10s %12s %12s\n', 'age', 'mass', 'mean mu', 'mean MRPK', 'sd log MRPK');
    for aa = 1:min(size(INV_AGE_COH,2), 20)
        if isfinite(R(e).age_mu(aa))
            fprintf('  %5d %11.5f %10.4f %12.4f %12.4f\n', ...
                aa, R(e).age_mass(aa), R(e).age_mu(aa), R(e).age_mrpk(aa), R(e).age_sdm(aa));
        end
    end
end

AK   = size(INV_AGE_COH, 2);
yidx = 1:min(5,AK);
oidx = min(15,AK):min(AK,40);
my = zeros(1,NE); mo = zeros(1,NE); sy = zeros(1,NE); so = zeros(1,NE);
for e = 1:NE
    my(e) = mean(R(e).age_mrpk(yidx),'omitnan');
    mo(e) = mean(R(e).age_mrpk(oidx),'omitnan');
    sy(e) = mean(R(e).age_sdm(yidx),'omitnan');
    so(e) = mean(R(e).age_sdm(oidx),'omitnan');
end

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SELECTION AND ALLOCATION -- CF friction against CF frictionless\n');
fprintf('================================================================\n\n');
fprintf('  %-42s %14s %14s\n', 'statistic', 'friction', 'frictionless');
fprintf('  %-42s %14.5f %14.5f\n', 'output Y', R(1).Y, R(2).Y);
fprintf('  %-42s %14.5f %14.5f\n', 'TFP', R(1).TFP, R(2).TFP);
fprintf('  %-42s %14.4f %14.4f\n', 'entrepreneur share', R(1).ent, R(2).ent);
fprintf('  %-42s %13.1f%% %13.1f%%\n', 'constrained entrepreneurs', 100*R(1).constr, 100*R(2).constr);
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), all firms', R(1).sd_logMRPK, R(2).sd_logMRPK);
fprintf('  %-42s %14.4f %14.4f\n', 'mean MRPK, ages 1-5', my(1), my(2));
fprintf('  %-42s %14.4f %14.4f\n', 'mean MRPK, ages 15+', mo(1), mo(2));
fprintf('  %-42s %+13.1f%% %+13.1f%%\n', 'young MRPK premium (Midrigan-Xu: +73%)', ...
    100*(my(1)/mo(1)-1), 100*(my(2)/mo(2)-1));
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), ages 1-5', sy(1), sy(2));
fprintf('  %-42s %14.4f %14.4f\n', 'sd(log MRPK), ages 15+', so(1), so(2));
fprintf('  %-42s %14.4f %14.4f\n', 'mean mu, ages 1-5', mean(R(1).age_mu(yidx),'omitnan'), mean(R(2).age_mu(yidx),'omitnan'));

fprintf('\n');
fprintf('================================================================\n');
fprintf('  MISALLOCATION LOSS  100*log(TFP frictionless / TFP friction)\n');
fprintf('  Midrigan convention. Both economies here lack used capital, so this is\n');
fprintf('  the loss the constraint causes WITHOUT a second-hand market to soften it.\n');
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

fprintf('\n  CHECK -- recomputed here at those same cleared prices:\n');
fprintf('    TFP %.6f and %.6f  ->  loss %.4f   (difference %.2e)\n', ...
    R(1).TFP, R(2).TFP, loss_self, loss_self - loss_saved);
if abs(loss_self - loss_saved) > 0.05
    fprintf('  ^ THESE DISAGREE. Trust the authoritative block and re-run nofriction.m.\n');
else
    fprintf('    Agreement confirms the analysis sits on the GE-cleared equilibrium.\n');
end

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  CHEN (2023) DIAGNOSTICS -- computed WITHIN this economy\n');
fprintf('================================================================\n');

AKc = size(INV_AGE_COH, 2);
mrpk_star = R(2).age_mrpk(min(15,AKc));
if ~isfinite(mrpk_star), mrpk_star = mean(R(2).age_mrpk,'omitnan'); end
tolsf = 0.05;
psf = find(R(1).age_mrpk <= mrpk_star*(1+tolsf) & isfinite(R(1).age_mrpk), 1, 'first');
fprintf('\n  PERIODS TO SELF-FINANCE\n');
fprintf('    unconstrained benchmark: MRPK = %.4f (this economy''s frictionless twin)\n', mrpk_star);
fprintf('    first age at which mean MRPK falls within %.0f%% of it: ', 100*tolsf);
if isempty(psf)
    fprintf('NOT within %d years\n', AKc);
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
fprintf('    %-10s %14s %14s %14s\n', '', 'K friction', 'K frictionless', 'relative');
for q = 1:4
    fprintf('    Q%-9d %14.5f %14.5f %+13.3f\n', q, R(1).Kq(q), R(2).Kq(q), ...
        R(2).Kq(q)/max(R(1).Kq(q),1e-12) - 1);
end
fprintf('    Q1 = lowest TFP.  Chen: -0.300, -0.269, -0.134, +0.216.\n');

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
fprintf('\n    Poor talented agents enter at %.1f%% of the rate perfect credit would give.\n', ...
    100*R(1).pt_rate/max(R(2).pt_rate,1e-12));

fprintf('\n  This file covers the CF alone. The benchmark-versus-CF comparison, and\n');
fprintf('  hence what used capital is worth, is produced by paper_tables.m.\n');

A = struct('econ',{R}, 'misalloc_loss',loss_saved, 'misalloc_self',loss_self, 'theta_cal',theta_cal);
save('analysis_cf_result.mat','A');

f = figure('Visible','off','Position',[100 100 900 380]);
nplot = min(size(INV_AGE_COH,2), 25);
subplot(1,2,1);
plot(1:nplot, R(1).age_mu(1:nplot), '-o', 'LineWidth', 1.6, 'MarkerSize', 4); hold on;
plot(1:nplot, R(2).age_mu(1:nplot), '-s', 'LineWidth', 1.6, 'MarkerSize', 4);
xlabel('firm age (years since entry)');
ylabel('shadow cost \mu');
legend({'friction','frictionless'}, 'Location','northeast');
title('CF: collateral shadow cost by age');
grid on;
subplot(1,2,2);
plot(1:nplot, R(1).age_sdm(1:nplot), '-o', 'LineWidth', 1.6, 'MarkerSize', 4); hold on;
plot(1:nplot, R(2).age_sdm(1:nplot), '-s', 'LineWidth', 1.6, 'MarkerSize', 4);
xlabel('firm age (years since entry)');
ylabel('sd(log MRPK)');
legend({'friction','frictionless'}, 'Location','northeast');
title('CF: MRPK dispersion by age');
grid on;
print(f, 'fig_analysis_cf.png', '-dpng', '-r150');
close(f);

fprintf('\n  [saved] analysis_cf_result.mat\n');
fprintf('  [saved] fig_analysis_cf.png\n');

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
