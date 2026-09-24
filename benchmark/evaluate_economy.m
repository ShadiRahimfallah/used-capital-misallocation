function y = evaluate_economy(x)
global rho alpha eta beta delta W r smin smax fspacee fspacew fspace ce cw ...
    egrid P smine sminw smaxw theta xai F cxe cxw cse csw Lbar phi gamma se ...
    Rnew eps_model productivity_weight Perg k exit_prob smaxe qu zeta kappa ...
    pu calib_mode print_x0_manual_style ygrid Py Ny rho_y sigma_y Ey piy ...
    wage_gini_target sigma_crra lambda_u sigma_perm eta_p omega ...
    fast_diagnostic xsim_prev_for_amax gaps_calibration_last WARM_CE WARM_CW ...
    VF_WARMSTART_XCALL cf_gamma_override cf_phi_override ...
    hybrid_use_nested_clearing hybrid_inner_clearing_skip ...
    hybrid_nested_price_lb hybrid_nested_price_ub ...
    ergodic_skip_legacy_moments KnKu_ratio psi_u cf_cf1_zero_used_prices ...
    cf_force_qu_eq_Rnew cf_no_used_dynamics GRID_LEGACY amax_calib_mult ...
    AMAX_OVERRIDE VF_NA VF_CURV VF_N_OUTER VF_N_PICARD VF_N_HOWARD DtoY ...
    beta_size_data MOMENTS_FROM_INVARIANT equilibrium_cf_eval ...
    DROP_CLEARING_WEIGHTS firm_exit_target top10_emp_share_target ...
    nested_clearing_tol nested_clearing_tol_Knew calib_moment_weights ...
    R_TARGET WGINI_USE_TRIM

if nargin < 1
    error(['objective_calibration requires input x = [W; r; Rnew; entry_cost_slot; theta; se; rho; kappa; pw; sigma_y; eta; exit_prob]. ' ...
        'Run start_calibration to launch calibration.']);
end
if isempty(x) || numel(x) < 11
    error('evaluate_economy requires x with >=11 (legacy) or 12 (+eta_span) elements. Got %d.', numel(x));
end

type     = [];
typeold  = [];
ent_mask = [];

ce      = [];
cxe     = [];
cse     = [];
cw      = [];
cxw     = [];
csw     = [];
fspacee = [];
fspacew = [];

if isempty(beta), beta = 0.865; end
alpha      = 2/3;

delta = 0.06;
Lbar  = 1;
xai   = 0;
if isempty(F), F = 0.00; end
eps_model = 1e-6;
spliorder = 1;

if isempty(sigma_crra), sigma_crra = 1.5; end

if isempty(rho_y), rho_y = 0.95; end
if isempty(Ny),    Ny    = 7;    end
if isempty(wage_gini_target), wage_gini_target = 0.372; end

gamma = 2;
if ~isempty(cf_gamma_override) && ~any(isnan(cf_gamma_override(:)))
    gamma = cf_gamma_override(1);
end

if ~isempty(cf_phi_override) && ~any(isnan(cf_phi_override(:)))
    phi = cf_phi_override(1);
elseif isempty(phi) || ~isscalar(phi) || ~isfinite(phi) || ~(phi > 0 && phi <= 1)
    phi = 0.75;
end

if ~exist('calib_mode','var') || isempty(calib_mode)
    calib_mode = true;
end

W                   = x(1);
r                   = x(2);
Rnew                = x(3);
entry_cost_slot     = 0;
theta               = x(5);
productivity_weight = 1.0;

se  = 0;
rho = 0;

if numel(x) >= 12
    eta       = max(min(x(6), 0.95), 0.55);
    eta_p     = x(7);
    omega     = x(8);
    kappa     = x(9);
    sigma_y   = x(10);
    exit_prob = x(11);
    lambda_u  = x(12);
else
    eta       = 0.79;
    eta_p     = x(6);
    omega     = x(7);
    kappa     = x(8);
    sigma_y   = x(9);
    exit_prob = x(10);
    lambda_u  = x(11);
end
sigma_perm = 0.0;

if isempty(hybrid_use_nested_clearing), hybrid_use_nested_clearing = false; end
if isempty(hybrid_inner_clearing_skip), hybrid_inner_clearing_skip = false; end
if hybrid_use_nested_clearing && ~hybrid_inner_clearing_skip && calib_mode && numel(x) >= 11
    if isempty(hybrid_nested_price_lb) || numel(hybrid_nested_price_lb) ~= 3
        hybrid_nested_price_lb = [0.20; 0.025; 0.095];
    end
    if isempty(hybrid_nested_price_ub) || numel(hybrid_nested_price_ub) ~= 3
        hybrid_nested_price_ub = [1.20; 0.070; 0.160];
    end
    if numel(x) >= 12
        struct_p = x(4:12);
    else
        struct_p = x(4:11);
    end
    try
        [W, r, Rnew, ~, ~] = clear_markets( ...
            struct_p, W, r, Rnew, hybrid_nested_price_lb(:), hybrid_nested_price_ub(:));
    catch
    end
end

if isempty(KnKu_ratio), KnKu_ratio = 2.5699; end
zeta      = (delta + kappa) / KnKu_ratio;
spread_rd = Rnew - r - delta;

if isempty(psi_u), psi_u = 0; end
Rused = r + delta + kappa * (1 - psi_u);

tol_cf1_pin = max(1e-9, 0.05 * max(abs(zeta), 1e-12));
cf1_qu_off = (~isempty(cf_cf1_zero_used_prices) && logical(cf_cf1_zero_used_prices)) ...
    || (lambda_u < 1e-8 && abs(spread_rd - zeta) <= tol_cf1_pin);

if cf1_qu_off
    qu = 0;
    pu = 0;
else
    qu = Rused * (1 - spread_rd / zeta);

    tol_z = 1e-10 * max(1, abs(zeta));
    if spread_rd > zeta + tol_z
        y                     = 100;
        gaps_calibration_last = 11 * ones(8, 1);
        fprintf(['PENALTY: implied qu<0 (Rnew-r-delta=%.6f > zeta=%.6f). ', ...
            'Used asset price P_used would be negative.\n'], spread_rd, zeta);
        return;
    end
    if qu <= 1e-6
        qu = 0;
        pu = 0;
    else
        pu = qu / Rused;
    end
end

if ~calib_mode
    fprintf('\n--- EVALUATING x ---\n');
    fprintf('  W=%.4f  r=%.4f  Rnew=%.4f  qu=%.4f  theta=%.4f\n', ...
        W, r, Rnew, qu, theta);
    fprintf('  eta_p=%.4f  omega=%.3f  pu=%.3f (=qu/Rused)  kappa=%.4f  zeta=%.4f  pw=%.3f  sigma_y=%.3f\n', ...
        eta_p, omega, pu, kappa, zeta, productivity_weight, sigma_y);
end

if Rnew <= r + delta + 1e-6
    y                     = 100;
    gaps_calibration_last = 11 * ones(8, 1);
    fprintf('PENALTY: Rnew=%.4f <= r+delta=%.4f\n', Rnew, r+delta);
    return;
end

qu_spread_penalty = 0;

if isempty(print_x0_manual_style), print_x0_manual_style = false; end
if print_x0_manual_style
    fprintf('\n[objective_calibration] x (each evaluation)\n');
    print_parameters(x);
end

if ~isempty(cf_force_qu_eq_Rnew) && cf_force_qu_eq_Rnew
    qu                = Rnew;
    pu                = 1.0;
    qu_spread_penalty = 0;
end

if ~isempty(cf_no_used_dynamics) && cf_no_used_dynamics
    qu                = 0;
    pu                = 0;
    qu_spread_penalty = 0;
end

ap_safe = max(eta_p, 1.10);
om_safe = max(min(omega, 0.999), 0.0);

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic, GRID_LEGACY);

assert(all(abs(sum(P,2) - 1) < 1e-10), '[egrid] P rows must sum to 1');
assert(abs(sum(Perg) - 1) < 1e-10,     '[egrid] Perg must sum to 1');
assert(abs(pi_z' * zP - 1) < 1e-8,     '[egrid] E[z]=1 normalisation failed');

egrid_original = egrid;

if exist('calib_mode','var') && ~calib_mode
    fprintf('  [egrid] Buera-Shin (2013) Pareto + omega-redraw: k=%d, eta_p=%.3f, omega=%.3f\n', ...
        k, eta_p, om_safe);
    fprintf('  [egrid] zP range:           min=%.3f  median=%.3f  max=%.3f\n', ...
        min(zP), median(zP), max(zP));
    fprintf('  [egrid] exp(egrid) range=%.0fx (max/min)\n', max(zP)/min(zP));
    fprintf('  [egrid] tail mass at top-2 nodes (F>0.999) = %.4f\n', sum(pi_z(end-1:end)));
end

[ygrid, Py, Ey, piy] = rouwenhorst_y(Ny, rho_y, sigma_y);

try
    solvee_static_min

    Amin = Aumin;
    amin = max(Aumin + exp(-7) , 0.001);

    amax = 80;
    if ~isempty(AMAX_OVERRIDE), amax = AMAX_OVERRIDE; end

    if amin >= amax - 0.01 || (amax - Aumin) < 0.05
        y                     = 100;
        gaps_calibration_last = 11 * ones(8, 1);
        fprintf('PENALTY: invalid grid amin=%.4f amax=%.4f (Aumin=%.4f)\n', amin, amax, Aumin);
        return;
    end

    smin  = [amin, 1, 1 ];
    smax  = [amax, k, Ny];
    smine = smin;
    smaxe = smax;

    if ~isempty(fast_diagnostic) && fast_diagnostic
        n    = [100, k, Ny];
        curv = .1;
    else
        n    = [100, k, Ny];
        curv = .1;
    end

    if ~isempty(VF_NA)   && VF_NA   >= 20, n(1) = round(VF_NA); end
    if ~isempty(VF_CURV) && VF_CURV >  0,  curv = VF_CURV;      end
    Amin_ref = Aumin;
    agrid    = nodeunif(n(1), (amin-Amin_ref).^curv, (amax-Amin_ref).^curv).^(1/curv) + Amin_ref;

    fspace = fundef({'spli', agrid,     0, spliorder}, ...
        {'spli', (1:1:k)',  0, 1}, ...
        {'spli', (1:1:Ny)', 0, 1});
    grid = funnode(fspace);

    if iscell(grid)
        s = gridmake(grid{:});
    else
        s = gridmake(grid);
    end
    ns = size(s, 1);

    if size(s, 2) ~= 3
        error(['evaluate_economy: state s has %d columns but fspace ' ...
            'is 3D (a, z, y). Check fspace setup. Sizes: agrid=%dx%d, k=%d, Ny=%d.'], ...
            size(s, 2), size(agrid,1), size(agrid,2), k, Ny);
    end

    solvee_static;

    fspace_E = fspace;
    grid_E   = funnode(fspace_E);
    if iscell(grid_E)
        s_E = gridmake(grid_E{:});
    else
        s_E = gridmake(grid_E);
    end
    ns_E = size(s_E, 1);

    Amin_w = 0;
    amin_w = max(Amin_w + exp(-15) + 0.00, 0.001);
    amax_w = 80;
    if ~isempty(AMAX_OVERRIDE), amax_w = AMAX_OVERRIDE; end
    smin   = [amin_w, 1, 1 ];
    smax   = [amax_w, k, Ny];
    sminw  = smin;
    smaxw  = smax;

    n_w     = [50, k, Ny]; curv_w = .1; Amin_ref_w = 0;
    agrid_w = nodeunif(n_w(1),(amin_w-Amin_ref_w).^curv_w,(amax_w-Amin_ref_w).^curv_w).^(1/curv_w)+Amin_ref_w;
    fspace_W = fundef({'spli', agrid_w,    0, spliorder}, ...
        {'spli', (1:1:k)',   0, 1}, ...
        {'spli', (1:1:Ny)',  0, 1});
    grid_W   = funnode(fspace_W);
    if iscell(grid_W)
        s_W = gridmake(grid_W{:});
    else
        s_W = gridmake(grid_W);
    end
    ns_W = size(s_W, 1);

    if size(s_W, 2) ~= 3
        error('evaluate_economy: worker state s_W has %d columns (expected 3)', size(s_W,2));
    end

    if ~isempty(fast_diagnostic) && fast_diagnostic
        n_outer  = 12;
        n_picard = 60;
    else
        n_outer  = 40;
        n_picard = 200;
    end

    if ~isempty(VF_N_OUTER),  n_outer  = VF_N_OUTER;  end
    if ~isempty(VF_N_PICARD), n_picard = VF_N_PICARD; end

    if isempty(VF_N_HOWARD), n_howard = 30; else, n_howard = VF_N_HOWARD; end
    fprintf('  [VF] solving (Na=%d, Nz=%d, Ny=%d => ns=%d per type)...\n', ...
        n(1), k, Ny, ns_E);

    if isempty(VF_WARMSTART_XCALL) || VF_WARMSTART_XCALL
        VFCACHE = fullfile(fileparts(mfilename('fullpath')), 'vf_warm_cache.mat');
        if isempty(WARM_CE) && isempty(WARM_CW) && exist(VFCACHE, 'file')
            try
                Wc = load(VFCACHE);
                if isfield(Wc,'WARM_CE') && isequal(size(Wc.WARM_CE), [ns_E, 2])
                    WARM_CE = Wc.WARM_CE;
                end
                if isfield(Wc,'WARM_CW') && isequal(size(Wc.WARM_CW), [ns_W, 2])
                    WARM_CW = Wc.WARM_CW;
                end
                if ~isempty(WARM_CE) || ~isempty(WARM_CW)
                    fprintf('  [VF] disk cache HIT (%s)\n', VFCACHE);
                end
            catch
                fprintf('  [VF] disk cache unreadable -- cold start.\n');
            end
        end
        if ~isempty(WARM_CE) && isequal(size(WARM_CE), [ns_E, 2])
            ce = WARM_CE;  fspacee = fspace_E;
        end
        if ~isempty(WARM_CW) && isequal(size(WARM_CW), [ns_W, 2])
            cw = WARM_CW;  fspacew = fspace_W;
        end
        if ~isempty(ce) || ~isempty(cw)
            fprintf('  [VF] warm start from previous call (ce %s, cw %s)\n', ...
                mat2str(~isempty(ce)), mat2str(~isempty(cw)));
        end
    end

    t_vf = tic;
    for outer = 1:n_outer
        fspace = fspace_E;

        if ~isempty(ce) && isequal(size(ce), [ns_E, 2])
            c = ce;
        else
            v = zeros(ns_E,1);
            c = funfitxy(fspace,s_E,[v,v]);
        end

        for i = 1:n_picard
            cnew = c;
            [v1,~,~,x_e] = saveBelmaxe(cnew,fspace,s_E);
            c(:,1) = v1;
            v2     = valfunc2e(c,fspace,s_E);
            c(:,2) = v2;
            for h = 1:n_howard
                chow   = c;
                v1h    = valfunc1e(c,fspace,s_E,x_e);
                c(:,1) = v1h;
                v2h    = valfunc2e(c,fspace,s_E);
                c(:,2) = v2h;
                if norm(c-chow)/norm(c) < 1e-6, break, end
            end
            if norm((c-cnew))/norm(c) < 1e-5, break, end
        end
        fprintf('  [VF outer %d] E Picard done in %d/%d iters (%.1fs)\n', ...
            outer, i, n_picard, toc(t_vf));

        c = vec(c);
        c = [c(1:ns_E), c(ns_E+1:2*ns_E)];
        [~,~,~,xe] = saveBelmaxe(c,fspace,s_E);
        cx      = funfitxy(fspace,s_E,xe);
        fspacee = fspace;
        ce      = c;
        cxe     = cx;
        [~,vee,vwe] = valfunc2e(c,fspace,s_E);
        cse = funfitxy(fspace,s_E,[vee,vwe]);

        if exist('cw','var') && ~isempty(cw)
            cw_old = cw;
        else
            cw_old = [];
        end

        fspace = fspace_W;

        if ~isempty(cw) && isequal(size(cw), [ns_W, 2])
            c = cw;
        else
            v = zeros(ns_W,1);
            c = funfitxy(fspace,s_W,[v,v]);
        end
        for i = 1:n_picard
            cnew = c;
            [v1,~,x_w] = saveBelmaxw(cnew,fspace,s_W);
            c(:,1) = v1;
            v2     = valfunc2w(c,fspace,s_W);
            c(:,2) = v2;
            for h = 1:n_howard
                chow   = c;
                v1h    = valfunc1w(c,fspace,s_W,x_w);
                c(:,1) = v1h;
                v2h    = valfunc2w(c,fspace,s_W);
                c(:,2) = v2h;
                if norm(c-chow)/norm(c) < 1e-6, break, end
            end
            if norm((c-cnew))/norm(c) < 1e-5, break, end
        end
        fprintf('  [VF outer %d] W Picard done in %d/%d iters (%.1fs)\n', ...
            outer, i, n_picard, toc(t_vf));

        c = vec(c);
        c = [c(1:ns_W), c(ns_W+1:2*ns_W)];
        [~,~,x_] = saveBelmaxw(c,fspace,s_W);
        cx      = funfitxy(fspace,s_W,x_);
        fspacew = fspace;
        cw      = c;
        cxw     = cx;
        [~,vew,vww] = valfunc2w(cw,fspacew,s_W);
        csw = funfitxy(fspacew,s_W,[vew,vww]);

        if ~isempty(cw_old) && all(size(cw_old) == size(cw))
            cross_change = norm(cw(:) - cw_old(:)) / max(norm(cw(:)), 1e-10);
            if cross_change < 1e-5
                break;
            end
        end
    end

    WARM_CE = ce;
    WARM_CW = cw;

    try
        VFCACHE = fullfile(fileparts(mfilename('fullpath')), 'vf_warm_cache.mat');
        save(VFCACHE, 'WARM_CE', 'WARM_CW');
    catch
        fprintf('  [VF] could not write the warm-start cache (continuing).\n');
    end

    fspace = fspacee;

    clear bel beljac c cx v1 v2 x_ cw_old cross_change

    if isempty(DtoY),      DtoY = 0.6000; end
    if isempty(KnKu_ratio), KnKu_ratio = 2.5699; end
    if isempty(beta_size_data), beta_size_data = -0.0269; end

    rng(42, 'twister');

    if isempty(MOMENTS_FROM_INVARIANT), MOMENTS_FROM_INVARIANT = true; end
    if MOMENTS_FROM_INVARIANT
        moments_from_invariant;

        if ~exist('beta_size_data','var')    || isempty(beta_size_data),     beta_size_data    = -0.0269; end
        if ~exist('qu_spread_penalty','var') || isempty(qu_spread_penalty),  qu_spread_penalty = 0;      end
        if isempty(equilibrium_cf_eval), equilibrium_cf_eval = false; end

        gap_Knew  = gap_Kn;
        gap_Labor = gap_L;
        gap_Kused = gap_Ku;
        gap_KnKu  = (KnKu_sim - KnKu_ratio) / max(KnKu_ratio, 1e-6);
        gap_DtoY  = (DtoY_sim - DtoY)        / max(DtoY,       1e-6);

        scale_beta_size     = 0.03;
        beta_size_model_sim = beta_size_model;
        gap_beta_size       = (beta_size_model_sim - beta_size_data) / scale_beta_size;

        gini_wage_score = gini_wage_model;
        if ~isempty(WGINI_USE_TRIM) && WGINI_USE_TRIM ...
                && exist('gini_wage_trim','var') && isfinite(gini_wage_trim)
            gini_wage_score = gini_wage_trim;
        end
        UNCOMPUTED_MOMENT_PENALTY = 3;

        if isfinite(gini_wage_score)
            gap_wageGini = (gini_wage_score - wage_gini_target) / max(wage_gini_target, 1e-6);
        else
            gap_wageGini = UNCOMPUTED_MOMENT_PENALTY;
            fprintf('  [WARN] wage Gini not computable -- scored as a %d-fold miss, not a match\n', UNCOMPUTED_MOMENT_PENALTY);
        end
        if isfinite(firm_exit_model)
            gap_exit = (firm_exit_model - firm_exit_target) / max(firm_exit_target, 1e-6);
        else
            gap_exit = UNCOMPUTED_MOMENT_PENALTY;
            fprintf('  [WARN] exit rate not computable -- scored as a %d-fold miss, not a match\n', UNCOMPUTED_MOMENT_PENALTY);
        end
        if isfinite(top10_share_model)
            gap_top10 = (top10_share_model - top10_emp_share_target) / max(top10_emp_share_target, 1e-6);
        else
            gap_top10 = UNCOMPUTED_MOMENT_PENALTY;
            fprintf('  [WARN] top10 share not computable -- scored as a %d-fold miss, not a match\n', UNCOMPUTED_MOMENT_PENALTY);
        end

        gaps_calibration_last = [gap_Knew; gap_Labor; gap_KnKu; gap_DtoY; ...
            gap_beta_size; gap_wageGini; gap_exit; ...
            gap_top10];
        y_vec         = gaps_calibration_last;
        current_share = ent_share_sim;

        if equilibrium_cf_eval
            clr_cap = max(-3, min(3, [gap_Knew; gap_Labor; gap_KnKu]));
            y       = sqrt(clr_cap(1)^2 + clr_cap(2)^2 + 0.25*clr_cap(3)^2);
            if abs(gap_Knew) > 0.35, y = y + 0.20 * (abs(gap_Knew) - 0.35); end
            y = y + qu_spread_penalty;
        else
            clr_cap = max(-3, min(3, [gap_Knew; gap_Labor; gap_KnKu]));
            tlp_h   = nested_clearing_tol;
            if isempty(tlp_h) || tlp_h <= 0, tlp_h = 0.01;
            end
            if ~isempty(nested_clearing_tol_Knew) && nested_clearing_tol_Knew > 0
                tol_kn_mkt = nested_clearing_tol_Knew;
            else
                tol_kn_mkt = 0;
            end
            gKn_mkt = sign(clr_cap(1)) * max(abs(clr_cap(1)) - tol_kn_mkt, 0);
            mkt_obj = sqrt(gKn_mkt^2 + clr_cap(2)^2 + 0.25*clr_cap(3)^2);

            if isempty(R_TARGET), R_TARGET = 0.031; end
            gap_r = (r - R_TARGET) / max(R_TARGET, 1e-6);

            mom_vec = [gap_DtoY; gap_beta_size; gap_wageGini; gap_exit; gap_top10];
            if ~isempty(calib_moment_weights) && numel(calib_moment_weights) == 6
                mom_vec = [mom_vec; gap_r];
                w_mom   = calib_moment_weights(:);
            elseif ~isempty(calib_moment_weights) && numel(calib_moment_weights) == 5
                w_mom = calib_moment_weights(:);
            else
                w_mom = [35; 20; 35; 30; 30];
            end
            mom_cap    = max(-3, min(3, mom_vec));
            moment_obj = sqrt(sum(mom_cap.^2 .* w_mom) / sum(w_mom));
            if ~isempty(DROP_CLEARING_WEIGHTS) && DROP_CLEARING_WEIGHTS
                y = moment_obj;
            else
                y = sqrt(mkt_obj^2 + moment_obj^2);
            end
            if abs(gap_Knew) > 0.35, y = y + 0.20 * (abs(gap_Knew) - 0.35); end
            if current_share < 0.02, y = y + 10.0; end
            if current_share > 0.85, y = y + 10.0; end
            y = y + qu_spread_penalty;
        end

        if ~isempty(calib_mode) && ~calib_mode
            fprintf('objective (TAN-FAST)          = %9.6f \n', y);
            fprintf('\n');
        end
        return;
    end

    error('evaluate_economy:legacyPath', ...
        ['The legacy panel-simulation path is not shipped. Set ' ...
         'MOMENTS_FROM_INVARIANT = true, which every driver in this ' ...
         'package already does.']);

    if ~exist('ent_share_sim','var'), ent_share_sim = 0.20; end

    if exist('xsim','var') && ~isempty(xsim)
        xsim_prev_for_amax = xsim;
        frac_at_top        = mean(xsim(:) > 0.99 * smaxe(1));
        if frac_at_top > 0.005
            fprintf('  [amax-WARN] %.2f%% of savings hit upper bound (smaxe=%.2f). ', ...
                100*frac_at_top, smaxe(1));
            fprintf('Wealth tail compressed; raise amax adaptively next iter.\n');
        end
    end

    if exist('Kn','var') && exist('Ku','var') && exist('ent_mask','var') && exist('A','var')
        mask_size   = ent_mask & (Kn + Ku > 1e-12);
        n_mask_size = sum(mask_size(:));
        if n_mask_size >= 5
            us_reg  = Ku(mask_size) ./ (Kn(mask_size) + Ku(mask_size));
            lnA_reg = log(A(mask_size));
            ok_size = isfinite(us_reg) & isfinite(lnA_reg);
            us_reg  = us_reg(ok_size);
            lnA_reg = lnA_reg(ok_size);
            if length(us_reg) >= 5
                p_size          = polyfit(lnA_reg, us_reg, 1);
                beta_size_model = p_size(1);
            else
                beta_size_model = 0;
            end
        else
            beta_size_model = 0;
        end

        if ~calib_mode
            fprintf('  beta_size=%.6f (N=%d), pu=%.4f (=qu/Rused, derived)\n', ...
                beta_size_model, n_mask_size, pu);
        end
    else
        beta_size_model = 0;
        fprintf('  WARNING: missing Kn/Ku — beta_size moment set to 0\n');
    end

    if exist('Kn','var') && exist('Ku','var') && exist('L','var') && exist('ent_mask','var') && exist('A','var')
        T_rec_sync   = size(ent_mask, 1);
        T_avg_sync   = min(30, T_rec_sync);
        t0_sync      = max(1, T_rec_sync - T_avg_sync + 1);
        ent_mask_win = ent_mask(t0_sync:end, :);
        Kn_w         = Kn(t0_sync:end, :);
        Ku_w         = Ku(t0_sync:end, :);
        L_w          = L(t0_sync:end, :);
        Debt_w       = Debt(t0_sync:end, :);
        Y_w          = Y(t0_sync:end, :);
        a_w          = a(t0_sync:end, :);

        n_obs            = numel(ent_mask_win);
        Knew_demand_sim  = sum(Kn_w(ent_mask_win(:))) / n_obs;
        Kused_demand_sim = sum(Ku_w(ent_mask_win(:))) / n_obs;
        L_demand_sim     = sum(L_w(ent_mask_win(:))) / n_obs;
        KnKu_sim         = sum(Kn_w(ent_mask_win(:))) / max(sum(Ku_w(ent_mask_win(:))), 1e-10);

        DtoY_sim     = sum(Debt_w(ent_mask_win(:))) / max(sum(Y_w(ent_mask_win(:))), 1e-10);
        DtoY_sim_net = NaN;
        if exist('Debt_net','var') && ~isempty(Debt_net)
            Debt_net_w   = Debt_net(t0_sync:end, :);
            DtoY_sim_net = sum(Debt_net_w(ent_mask_win(:))) / max(sum(Y_w(ent_mask_win(:))), 1e-10);
        end

        Dtot_sim         = sum(a_w(:)) / n_obs;

        Knew_supply_sim  = Dtot_sim / (1 + zeta*qu / (Rused*(delta+kappa)));
        Kused_supply_sim = zeta * Knew_supply_sim / (delta + kappa);

        if exist('y_idx','var') && ~isempty(ygrid)
            y_w            = y_idx(t0_sync:end, :);
            worker_mask    = ~ent_mask_win(:);
            L_supply_sim   = sum(ygrid(y_w(worker_mask))) / n_obs;
            L_supply_heads = sum(worker_mask) / n_obs;
        else
            L_supply_sim   = sum(~ent_mask_win(:)) / n_obs;
            L_supply_heads = L_supply_sim;
        end
    else
        Knew_demand_sim  = 0;
        Kused_demand_sim = 0;
        L_demand_sim     = 0;
        Knew_supply_sim  = 1;
        Kused_supply_sim = 1;
        L_supply_sim     = 1;
        KnKu_sim         = 2.95;
        DtoY_sim         = 0.64;
    end

    acorY_sim           = acorY;
    ent_share_sim_save  = ent_share_sim;
    beta_size_model_sim = beta_size_model;

    L_demand_sim_pre     = L_demand_sim;
    L_supply_sim_pre     = L_supply_sim;
    Knew_demand_sim_pre  = Knew_demand_sim;
    Knew_supply_sim_pre  = Knew_supply_sim;
    Kused_demand_sim_pre = Kused_demand_sim;
    Kused_supply_sim_pre = Kused_supply_sim;

    ergodic_histogram_from_sim_panel = true;
    ergodic_mode                     = true;
    erg_skip_save                    = [];
    if exist('ergodic_skip_legacy_moments', 'var')
        erg_skip_save = ergodic_skip_legacy_moments;
    end
    ergodic_skip_legacy_moments = true;
    error('evaluate_economy:legacyPath', ...
        ['The legacy panel-simulation path is not shipped (ergodic.m is ' ...
         'absent). Set MOMENTS_FROM_INVARIANT = true, which every driver ' ...
         'in this package already does.']);
    if isempty(erg_skip_save)
        clear ergodic_skip_legacy_moments
    else
        ergodic_skip_legacy_moments = erg_skip_save;
    end

    gap_Knew_erg  = gap_Kn;
    gap_Kused_erg = gap_Ku;
    gap_Labor_erg = gap_L;

    Knew_dem_erg  = Knew_demand_sim;
    Knew_sup_erg  = Knew_supply_sim;
    Kused_dem_erg = Kused_demand_sim;
    Kused_sup_erg = Kused_supply_sim;
    L_dem_erg     = L_demand_sim;
    L_sup_erg     = L_supply_sim;
    KnKu_erg      = KnKu_sim;
    DtoY_erg      = DtoY_sim;

    acorY         = acorY_sim;
    ent_share_sim = ent_share_sim_save;

    ent_sim_print = ent_share_sim;
    ent_erg_print = NaN;
    if exist('ne_erg','var') && ~isempty(ne_erg), ent_erg_print = sum(ne_erg); end

    gap_L_sim  = (L_demand_sim_pre     - L_supply_sim_pre)     / max(L_supply_sim_pre,     1e-10);
    gap_Kn_sim = (Knew_demand_sim_pre  - Knew_supply_sim_pre)  / max(Knew_supply_sim_pre,  1e-10);
    gap_Ku_sim = (Kused_demand_sim_pre - Kused_supply_sim_pre) / max(Kused_supply_sim_pre, 1e-10);

    fprintf('  ====================  SIMULATE vs ERGODIC  ====================\n');
    fprintf('                          simulate    ergodic     diff\n');
    fprintf('  ent_share         :  %10.4f  %10.4f  %+9.4f\n', ...
        ent_sim_print,    ent_erg_print,    ent_sim_print - ent_erg_print);
    fprintf('  L_demand          :  %10.4f  %10.4f  %+9.4f\n', ...
        L_demand_sim_pre, L_dem_erg,        L_demand_sim_pre - L_dem_erg);
    fprintf('  L_supply          :  %10.4f  %10.4f  %+9.4f\n', ...
        L_supply_sim_pre, L_sup_erg,        L_supply_sim_pre - L_sup_erg);
    fprintf('  gap_L (%%)         :  %+9.2f%%  %+9.2f%%  %+8.2fpp\n', ...
        100*gap_L_sim, 100*gap_Labor_erg, 100*(gap_L_sim - gap_Labor_erg));
    fprintf('  gap_Knew (%%)      :  %+9.2f%%  %+9.2f%%  %+8.2fpp\n', ...
        100*gap_Kn_sim, 100*gap_Knew_erg, 100*(gap_Kn_sim - gap_Knew_erg));
    fprintf('  ===============================================================\n');

    if ~calib_mode
        fprintf('\n  --- HYBRID: clearing gaps from ergodic ---\n');
        fprintf('  gap_Knew  (sim=%+.1f%%, erg=%+.1f%%)\n', 100*gap_Kn_sim, 100*gap_Knew_erg);
        fprintf('  gap_Labor (sim=%+.1f%%, erg=%+.1f%%)\n', 100*gap_L_sim,  100*gap_Labor_erg);
    end

    gap_Knew  = gap_Knew_erg;
    gap_Kused = gap_Kused_erg;
    gap_Labor = gap_Labor_erg;
    gap_KnKu  = (KnKu_erg - KnKu_ratio) / max(KnKu_ratio, 1e-6);
    gap_DtoY  = (DtoY_sim - DtoY) / max(DtoY, 1e-6);

    scale_beta_size = 0.03;
    gap_beta_size   = (beta_size_model_sim - beta_size_data) / scale_beta_size;

    gini_wage_model = NaN;
    if exist('y_idx','var') && ~isempty(ygrid) && ~isempty(Ey) && Ey > 0 ...
            && exist('ent_mask','var')
        if ndims(ent_mask) == 2 && size(ent_mask,1) > 1 && size(ent_mask,2) > 1
            T_panel  = size(ent_mask, 1);
            Nf_panel = size(ent_mask, 2);
            panel_2d = true;
        else
            T_panel  = NaN;
            Nf_panel = NaN;
            panel_2d = false;
        end

        if panel_2d
            gini_per_t = NaN(T_panel, 1);
            for tt = 1:T_panel
                wm_t = ~ent_mask(tt, :);
                if sum(wm_t) > 50
                    y_workers_t = ygrid(y_idx(tt, wm_t));
                    y_workers_t = y_workers_t(y_workers_t > 0);
                    if numel(y_workers_t) >= 50
                        qlo       = quantile(y_workers_t, 0.01);
                        qhi       = quantile(y_workers_t, 0.99);
                        y_trimmed = y_workers_t(y_workers_t >= qlo & y_workers_t <= qhi);
                        if numel(y_trimmed) >= 30
                            gini_per_t(tt) = local_gini(y_trimmed);
                        end
                    end
                end
            end
            gini_per_t = gini_per_t(isfinite(gini_per_t));
            if ~isempty(gini_per_t)
                gini_wage_model = mean(gini_per_t);
            end
        else
            worker_mask  = ~ent_mask(:);
            if sum(worker_mask) >= 50
                wages = ygrid(y_idx(worker_mask));
                wages = wages(wages > 0);
                if numel(wages) >= 50
                    qlo   = quantile(wages, 0.01);
                    qhi   = quantile(wages, 0.99);
                    wages = wages(wages >= qlo & wages <= qhi);
                    if numel(wages) >= 30
                        gini_wage_model = local_gini(wages);
                    end
                end
            end
        end
    end
    if isfinite(gini_wage_model)
        gap_wageGini = (gini_wage_model - wage_gini_target) / max(wage_gini_target, 1e-6);
    else
        gap_wageGini = 0;
    end

    if isempty(firm_exit_target),  firm_exit_target  = 0.09782; end

    if exist('type','var') && ~isempty(type) && size(type, 1) >= 2
        type_prev   = type(1:end-1, :);
        type_curr   = type(2:end,   :);
        exits_total = sum(type_prev(:) == 1 & type_curr(:) == 0);
        base_total  = sum(type_prev(:) == 1);
        if base_total > 0
            firm_exit_model = exits_total / base_total;
        else
            firm_exit_model = NaN;
        end
    else
        firm_exit_model = NaN;
    end
    if isfinite(firm_exit_model)
        gap_exit = (firm_exit_model - firm_exit_target) / max(firm_exit_target, 1e-6);
    else
        gap_exit = 0;
    end

    if exist('L','var') && exist('ent_mask','var') && sum(ent_mask(:)) > 10
        firm_size_model = mean(L(ent_mask(:)));
    else
        firm_size_model = NaN;
    end

    if isempty(top10_emp_share_target), top10_emp_share_target = 0.63210; end

    if exist('L','var') && exist('ent_mask','var') && sum(ent_mask(:)) > 100
        L_ent_vec = L(ent_mask(:));
        L_ent_vec = L_ent_vec(L_ent_vec > 0);
        if numel(L_ent_vec) > 10
            L_sorted          = sort(L_ent_vec, 'descend');
            n_top10           = max(1, ceil(numel(L_sorted) * 0.10));
            top10_share_model = sum(L_sorted(1:n_top10)) / sum(L_sorted);
            gap_top10         = (top10_share_model - top10_emp_share_target) / max(top10_emp_share_target, 1e-6);
        else
            top10_share_model = NaN;
            gap_top10         = 0;
        end
    else
        top10_share_model = NaN;
        gap_top10         = 0;
    end

    gaps_calibration_last = [gap_Knew; gap_Labor; gap_KnKu; gap_DtoY; ...
        gap_beta_size; gap_wageGini; gap_exit; ...
        gap_top10];

    current_share = ent_share_sim;

    y_vec = [gap_Knew; gap_Labor; gap_KnKu; gap_DtoY; ...
        gap_beta_size; gap_wageGini; gap_exit; ...
        gap_top10];

    if isempty(equilibrium_cf_eval), equilibrium_cf_eval = false; end

    if equilibrium_cf_eval
        clr_cap = max(-3, min(3, [gap_Knew; gap_Labor; gap_KnKu]));
        y       = sqrt(clr_cap(1)^2 + clr_cap(2)^2 + 0.25*clr_cap(3)^2);
        if abs(gap_Knew) > 0.35
            y = y + 0.20 * (abs(gap_Knew) - 0.35);
        end
        y = y + qu_spread_penalty;
    else
        clr_cap = max(-3, min(3, [gap_Knew; gap_Labor; gap_KnKu]));
        tlp_h   = nested_clearing_tol;
        if isempty(tlp_h) || tlp_h <= 0, tlp_h = 0.01;
        end
        if ~isempty(nested_clearing_tol_Knew) && nested_clearing_tol_Knew > 0
            tol_kn_mkt = nested_clearing_tol_Knew;
        else
            tol_kn_mkt = 0;
        end
        gKn_mkt = sign(clr_cap(1)) * max(abs(clr_cap(1)) - tol_kn_mkt, 0);
        mkt_obj = sqrt(gKn_mkt^2 + clr_cap(2)^2 + 0.25*clr_cap(3)^2);
        mom_vec = [gap_DtoY; gap_beta_size; gap_wageGini; gap_exit; gap_top10];
        if ~isempty(calib_moment_weights) && numel(calib_moment_weights) == 5
            w_mom = calib_moment_weights(:);
        else
            w_mom = [35; 20; 35; 30; 30];
        end
        mom_cap    = max(-3, min(3, mom_vec));
        moment_obj = sqrt(sum(mom_cap.^2 .* w_mom) / sum(w_mom));
        if ~isempty(DROP_CLEARING_WEIGHTS) && DROP_CLEARING_WEIGHTS
            y = moment_obj;
        else
            y = sqrt(mkt_obj^2 + moment_obj^2);
        end
        if abs(gap_Knew) > 0.35
            y = y + 0.20 * (abs(gap_Knew) - 0.35);
        end

        if current_share < 0.02
            y = y + 10.0;
            fprintf('[DEGENERATE] ent_share=%.4f < 0.02 — entry shut down. Penalty +10.\n', current_share);
        end
        if current_share > 0.85
            y = y + 10.0;
            fprintf('[DEGENERATE] ent_share=%.4f > 0.85 — everyone is entrepreneur. Penalty +10.\n', current_share);
        end
        if exist('L_supply_sim_pre','var') && L_supply_sim_pre < 0.05
            fprintf('[DEGENERATE] L_supply=%.4f < 0.05 — almost no workers in pool.\n', L_supply_sim_pre);
        end
        if exist('L_demand_sim_pre','var') && L_demand_sim_pre < 0.01
            fprintf('[DEGENERATE] L_demand=%.4f < 0.01 — firms hire almost no labor.\n', L_demand_sim_pre);
        end

        y = y + qu_spread_penalty;
    end

    if calib_mode
        if equilibrium_cf_eval
            error(['evaluate_economy: calib_mode=true with equilibrium_cf_eval ', ...
                '(equilibrium_cf) is unsupported; nested clearing must call this file directly.']);
        end

        sum_w = sum(w_mom);

        mom_contrib = mom_cap.^2 .* w_mom;
        if exist('gini_wage_score','var') && isfinite(gini_wage_score), wg_rep = gini_wage_score; else, wg_rep = gini_wage_model; end

        fprintf('\n');
        fprintf('PARAMETERS\n');
        fprintf('\n');
        fprintf('W                             = %9.4f \n', W);
        fprintf('r                             = %9.4f \n', r);
        fprintf('Rnew                          = %9.4f \n', Rnew);
        fprintf('theta                         = %9.4f \n', theta);
        fprintf('eta_span                      = %9.4f \n', eta);
        fprintf('eta_p                         = %9.4f \n', eta_p);
        fprintf('omega                         = %9.4f \n', omega);
        fprintf('kappa                         = %9.5f \n', kappa);
        fprintf('sigma_y                       = %9.4f \n', sigma_y);
        fprintf('exit_prob                     = %9.4f \n', exit_prob);
        fprintf('lambda_u                      = %9.4f \n', lambda_u);

        fprintf('\n');
        fprintf('MARKET CLEARING\n');
        fprintf('\n');
        disp('                                1st col: demand,   2nd col: supply,   3rd col: gap %');
        fprintf('\n');
        fprintf('Knew                          = %9.5f  %9.5f  %+9.1f \n', Knew_dem_erg, Knew_sup_erg, 100*gap_Knew_erg);
        fprintf('Labor                         = %9.5f  %9.5f  %+9.1f \n', L_dem_erg,    L_sup_erg,    100*gap_Labor_erg);
        fprintf('KnKu (model/target) [CLEARING]= %9.4f  %9.4f  %+9.1f \n', KnKu_erg,     KnKu_ratio,   100*gap_KnKu);

        fprintf('\n');
        fprintf('MOMENTS -- 5 TARGETED\n');
        fprintf('\n');
        disp('                                1st col: model,    2nd col: data,    3rd col: gap %,   4th col: cost');
        fprintf('\n');
        fprintf('DtoY (gross)                  = %9.4f  %9.4f  %+9.1f  %9.4f \n', ...
            DtoY_sim, DtoY, 100*gap_DtoY, mom_contrib(1)/sum_w);
        fprintf('beta_size (3rd col = scaled)  = %9.4f  %9.4f  %+9.4f  %9.4f \n', ...
            beta_size_model_sim, beta_size_data, gap_beta_size, mom_contrib(2)/sum_w);
        fprintf('wage Gini (scored)            = %9.4f  %9.4f  %+9.1f  %9.4f \n', ...
            wg_rep, wage_gini_target, 100*gap_wageGini, mom_contrib(3)/sum_w);
        fprintf('exit rate                     = %9.4f  %9.4f  %+9.1f  %9.4f \n', ...
            firm_exit_model, firm_exit_target, 100*gap_exit, mom_contrib(4)/sum_w);
        fprintf('top10 employment share        = %9.4f  %9.4f  %+9.1f  %9.4f \n', ...
            top10_share_model, top10_emp_share_target, 100*gap_top10, mom_contrib(5)/sum_w);

        fprintf('\n');
        fprintf('DIAGNOSTICS (untargeted)\n');
        fprintf('\n');
        fprintf('entrepreneur share            = %9.4f \n', ent_share_sim);
        fprintf('firm size (mean L)            = %9.4f \n', firm_size_model);
        fprintf('P_used = qu/Rused (derived)   = %9.4f \n', pu);
        if exist('DtoY_sim_net','var') && ~isnan(DtoY_sim_net)
            fprintf('DtoY (net, savers offset)     = %9.4f \n', DtoY_sim_net);
        end
        fprintf('wage Gini (untrimmed)         = %9.4f \n', gini_wage_model);

        fprintf('\n');
        fprintf('MOMENT WEIGHT TABLE (5 moments)\n');
        fprintf('\n');
        disp('                                1st col: weight,   2nd col: |gap| %,  3rd col: weighted cost');
        fprintf('\n');
        moment_lbls = {'DtoY','beta_size','wage Gini','exit rate','top10 emp share'};
        for ii = 1:5
            fprintf('%-30s= %9d  %9.1f  %9.4f \n', ...
                moment_lbls{ii}, w_mom(ii), 100*abs(mom_cap(ii)), mom_contrib(ii)/sum_w);
        end

        fprintf('\n');
        fprintf('OBJECTIVE\n');
        fprintf('\n');
        fprintf('mkt_obj (clearing)            = %9.4f \n', mkt_obj);
        fprintf('moment_obj (5 moments)        = %9.4f \n', moment_obj);
        fprintf('obj (total)                   = %9.4f \n', y);
        fprintf('\n');
    else
        fprintf('\n');
        fprintf('=== ALL GAPS (hybrid: clearing=ergodic, moments=simulate) ===\n');
        fprintf('\n');
        fprintf('CLEARING GAPS [ergodic]\n');
        fprintf('\n');
        fprintf('gap Knew                      = %+8.2f%% \n', 100*gap_Knew);
        fprintf('gap Labor                     = %+8.2f%% \n', 100*gap_Labor);
        fprintf('\n');
        fprintf('MOMENTS [simulate]\n');
        fprintf('\n');
        disp('                                1st col: model,    2nd col: data,    3rd col: gap %');
        fprintf('\n');
        fprintf('KnKu                          = %9.4f  %9.4f  %+9.2f \n', KnKu_sim, KnKu_ratio, 100*gap_KnKu);
        fprintf('DtoY                          = %9.4f  %9.4f  %+9.2f \n', DtoY_sim, DtoY,       100*gap_DtoY);
        fprintf('beta_size (3rd col = level)   = %9.4f  %9.4f  %+9.4f \n', ...
            beta_size_model_sim, beta_size_data, beta_size_model_sim - beta_size_data);

        fprintf('\n');
        fprintf('DERIVED\n');
        fprintf('\n');
        fprintf('zeta                          = %9.4f \n', zeta);
        fprintf('P_used = qu/Rused             = %9.4f \n', pu);

        fprintf('\n');
        fprintf('objective                     = %9.6f \n', y);
        if equilibrium_cf_eval
            fprintf('\n');
            fprintf(['  [equilibrium_cf] Scalar objective is clearing RMS only; ', ...
                'moment gaps vs Vietnam targets are not scored for CF.\n']);
        end
        fprintf('\n');
        fprintf('*******************************************************\n');
    end
catch ME
    y                     = 100;
    gaps_calibration_last = 11 * ones(8, 1);
    fprintf('PENALTY: %s\n', ME.message);

    fprintf('--- STACK TRACE ---\n');
    for k_st = 1:numel(ME.stack)
        fprintf('  [%d] %s  (line %d)  in %s\n', k_st, ...
            ME.stack(k_st).name, ME.stack(k_st).line, ME.stack(k_st).file);
    end
    fprintf('--- END STACK ---\n');
    return;
end

if isempty(equilibrium_cf_eval), equilibrium_cf_eval = false; end
if equilibrium_cf_eval
    return;
end
try
    if exist('best_equilibrium.mat','file')
        ld        = load('best_equilibrium.mat');
        prev_best = ld.best_obj;
    else
        prev_best = inf;
    end

    if y < prev_best
        best_save = struct();
        best_save.best_x = [W; r; Rnew];
        best_save.best_obj = y;
        best_save.best_iter = NaN;
        best_save.best_ent_share = current_share;
        best_save.theta = theta;
        best_save.eta = eta;
        best_save.kappa = kappa;
        best_save.sigma_y = sigma_y;
        best_save.exit_prob = exit_prob;
        best_save.zeta = zeta;
        best_save.pu = pu;
        best_save.ce = ce;
        best_save.cxe = cxe;
        best_save.cse = cse;
        best_save.cw = cw;
        best_save.cxw = cxw;
        best_save.csw = csw;
        best_save.fspacee = fspacee;
        best_save.fspacew = fspacew;
        best_save.smine = smine;
        best_save.smaxe = smaxe;
        best_save.sminw = sminw;
        best_save.smaxw = smaxw;
        save('best_equilibrium.mat', '-struct', 'best_save');
    end
catch
end
