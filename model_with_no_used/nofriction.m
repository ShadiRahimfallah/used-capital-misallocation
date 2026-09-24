clear;
clear global;
clc
format short

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
    psi_u WGINI_USE_TRIM INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT

FRICT_THETAS = [100000];

SEED_R    = 0.1015;
SEED_W_LO = 0.38;
SEED_W_HI = 0.80;
SEED_TOL  = 0.05;
SEED_MAX  = 7;

FRICT_SEED_W = 0.673115;
FRICT_SEED_R = 0.097422;

PROBE_W = [0.40 0.48 0.56 0.66 0.80];
PROBE_R = [0.060 0.100 0.140];

BIS_TOL       = 0.005;
BIS_MAX_EVALS = 160;

BASE_TOL = 0.015;

USE_BROYDEN = true;
SEED_J      = [ -9.0,  -6.1 ;
    -3.5, -22.6 ];
BRO_BOX       = [0.35 0.80 0.020 0.154];
BRO_MAX_EVALS = 30;

CKPT     = 'nofriction_frict_checkpoint.mat';
BASEFILE = 'nofriction_baseline.mat';

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
code_dir  = fileparts(script_dir);
repo_root = fileparts(code_dir);

if exist(fullfile(repo_root, 'CompEcon', 'CEtools'), 'dir')
    addpath(fullfile(repo_root, 'CompEcon', 'CEtools'), '-end');
elseif exist(fullfile(code_dir, 'CompEcon', 'CEtools'), 'dir')
    addpath(fullfile(code_dir, 'CompEcon', 'CEtools'), '-end');
end
if exist('menufun', 'file') ~= 2
    error('nofriction: menufun.m not on path.');
end

beta       = 0.853846484574;
alpha      = 2/3;
delta      = 0.06;
Lbar       = 1;
xai        = 0;
eps_model  = 1e-6;
sigma_crra = 1.5;

gamma = 2;
rho_y = 0.95;
Ny    = 7;

eta_span  = 0.79;
eta       = eta_span;
exit_prob = 0.0;

phi      = 1.0;
kappa    = 0.0;
zeta     = 0.0;
lambda_u = 0.0;
psi_u    = 0.0;
qu       = 0.0;
pu       = 0.0;

sigma_perm          = 0.0;
productivity_weight = 1.0;
se                  = 0;
rho                 = 0;

fast_diagnostic    = false;
LOCKOUT_EXOG_DEATH = true;
calib_mode         = true;
use_calib_sim_size = true;
amax_calib_mult    = 1.9;

hybrid_use_nested_clearing = false;
hybrid_inner_clearing_skip = false;

USE_INVARIANT_SEED = true;
SIM_NF             = 12000;
SIM_TINIT          = 2;

INV_NA                 = 3000;
INV_CURV               = 2.5;
INV_USE_EIGS           = true;
MOMENTS_FROM_INVARIANT = true;

nested_clearing_quiet    = false;
nested_clearing_max_iter = 60;
nested_clearing_tol      = 0.005;
nested_clearing_tol_Knew = 0.005;
calib_moment_weights     = [50; 35; 0; 30];
WGINI_USE_TRIM           = false;

hybrid_nested_price_lb = [0.35; 0.005; 0.04];
hybrid_nested_price_ub = [2.00; 0.99*(1/beta - 1); 0.26];

r_ceiling = 1/beta - 1;

DtoY                   = 0.6000;
wage_gini_target       = 0.372;
firm_exit_target       = 0.09782;
top10_emp_share_target = 0.63210;

parameters_cf_calib;
theta_cal = theta;
W_cal     = W;
r_cal     = r;
r_ceiling = 1/beta - 1;

PSTAMP = [theta_cal, eta_p, omega, sigma_y, beta, rho_y, F];

fprintf('  [nofriction] CF calibrated point read from parameters_cf_calib.m\n');
fprintf('    theta=%.6f eta_p=%.6f omega=%.6f sigma_y=%.6f beta=%.6f\n', ...
    theta_cal, eta_p, omega, sigma_y, beta);
fprintf('    base prices W=%.6f r=%.6f\n', W_cal, r_cal);

ap_safe = max(eta_p, 1.10);
om_safe = max(min(omega, 0.999), 0.0);

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, false);

[ygrid, Py, Ey, piy] = rouwenhorst_y(Ny, rho_y, sigma_y);

fprintf('================================================================\n');
fprintf('  NO-USED-CAPITAL CF -- ELIMINATE THE COLLATERAL CONSTRAINT\n');
fprintf('================================================================\n');
fprintf('  locked point: theta_cal=%.4f eta_p=%.4f omega=%.4f sigma_y=%.4f F=%.4f  (kappa/lambda_u removed)\n', ...
    theta_cal, eta_p, omega, sigma_y, F);
fprintf('  NOT recalibrated. Only prices are re-cleared. Rnew = r + delta throughout.\n');
fprintf('  frictionless thetas : %s   (matches the benchmark run)\n', mat2str(FRICT_THETAS));
fprintf('  r wall              : 1/beta-1 = %.4f\n', r_ceiling);
fprintf('================================================================\n');

base = [];
if exist(BASEFILE, 'file')
    bb = load(BASEFILE);
    if isfield(bb,'PSTAMP') && isfield(bb,'base') && numel(bb.PSTAMP)==numel(PSTAMP) ...
            && max(abs(bb.PSTAMP(:)'-PSTAMP)) < 1e-10
        base = bb.base;
        fprintf('\n  [baseline] loaded from %s (theta=%.4f, parameter stamp matches)\n', BASEFILE, base.theta);
    else
        fprintf('\n  [baseline] %s IGNORED -- built at different parameters. Re-solving.\n', BASEFILE);
    end
end
if isempty(base)
    fprintf('\n  [baseline] clearing theta=%.4f once and caching it (seed W=%.4f r=%.4f)...\n', ...
        theta_cal, W_cal, r_cal);
    theta  = theta_cal;
    W      = W_cal;
    r      = r_cal;
    Rnew   = r + delta;
    x_full = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    for wpass = 1:3
        evaluate_economy(x_full);
    end

    [Wc, rc, Rnewc, cvg, gb_out, nev_base] = clear_frictionless_broyden( ...
        theta, W, r, 'fd', BASE_TOL, BRO_MAX_EVALS, BRO_BOX);
    fprintf('  [baseline] Broyden used %d evals (tol %.1f%% -- see the note above).\n', ...
        nev_base, 100*BASE_TOL);

    x_full = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    evaluate_economy(x_full);
    gb   = gaps_calibration_last;
    base = aggregates();
    base.converged = cvg;
    base.gap_Knew  = 100*gb(1);
    base.gap_Labor = 100*gb(2);
    base.gap_KnKu  = 100*gb(3);
    base.LP        = base.Y / max(base.L, 1e-12);
    if ~cvg
        error(['nofriction: the CF BASELINE clear did not converge. Fix that first -- ' ...
            'everything below is measured against it.']);
    end
    if base.ent < 0.02 || base.ent > 0.90
        error(['nofriction: the CF baseline landed on the COLLAPSED occupational branch ' ...
            '(ent = %.3f). That is the known W knife-edge. Re-seed W and retry.'], base.ent);
    end
    save(BASEFILE, 'base', 'PSTAMP');
    fprintf('  [baseline] cleared and saved to %s\n', BASEFILE);
end

fprintf('  [baseline] theta=%.4f  W=%.4f r=%.4f Rnew=%.4f\n', base.theta, base.W, base.r, base.Rnew);
fprintf('             Y=%.5f  TFP=%.5f  ent=%.4f  constrained ents=%.1f%%  mean k/k*=%.4f\n', ...
    base.Y, base.TFP, base.ent, 100*base.constr, base.kk_mean);
fprintf('             residual gaps: Labor %+.2f%%  Knew %+.2f%%\n', base.gap_Labor, base.gap_Knew);
if max(abs([base.gap_Labor, base.gap_Knew])) > 0.5
    fprintf(['             ^ NOT cleared to 0.5%%. This economy''s capital market clears at the\n' ...
        '               EDGE OF AN r-CLIFF (dr = -0.0009 collapses ent 0.190 -> 0.009), so no\n' ...
        '               point on the healthy branch meets 0.5%%. Gains below are within-economy\n' ...
        '               RATIOS, so a ~1%% error in this baseline level moves the GAIN by ~1pp --\n' ...
        '               far less than the effect measured. Report this in the paper.\n']);
end

fr   = struct([]);
j0   = 1;
prev = [];
if exist(CKPT, 'file')
    ck = load(CKPT);
    if ~isfield(ck,'PSTAMP') || numel(ck.PSTAMP)~=numel(PSTAMP) ...
            || max(abs(ck.PSTAMP(:)'-PSTAMP)) >= 1e-10
        fprintf('\n  [resume] %s IGNORED -- different parameters or no stamp. Re-solving.\n', CKPT);
        ck = struct();
    end

    if isfield(ck,'fr') && isfield(ck,'FRICT_THETAS')
        nd = numel(ck.FRICT_THETAS);
        if nd <= numel(FRICT_THETAS) && isequal(ck.FRICT_THETAS(:)', FRICT_THETAS(1:nd))
            fr = ck.fr;
            j0 = numel(fr) + 1;
            ip = find([fr.converged], 1, 'last');
            if ~isempty(ip), prev = fr(ip); end
            fprintf('\n  [resume] %d frictionless theta(s) already done; continuing.\n', numel(fr));
        else
            fprintf('\n  [resume] checkpoint thetas differ from FRICT_THETAS -- starting over.\n');
        end
    end
end

for j = j0:numel(FRICT_THETAS)
    theta = FRICT_THETAS(j);
    fprintf('\n================================================================\n');
    fprintf('  FRICTIONLESS CLEAR:  theta = %.4g\n', theta);
    fprintf('================================================================\n');

    probe    = [];
    wlo      = [];
    whi      = [];
    rlo      = [];
    rhi      = [];
    nev_seed = 0;

    if isempty(prev) && ~isempty(FRICT_SEED_W)
        sW = FRICT_SEED_W;
        sR = FRICT_SEED_R;
        fprintf('\n  [seed] explicit frictionless seed W=%.6f r=%.6f -- wage bisection SKIPPED.\n', sW, sR);
        fprintf('         (set FRICT_SEED_W = [] at the top to bisect from scratch instead)\n');
    elseif isempty(prev)
        [sW, ~, ~, nev_seed] = seed_wage_bisect( ...
            theta, SEED_R, SEED_W_LO, SEED_W_HI, SEED_TOL, SEED_MAX);
        sR = SEED_R;
    else
        sW = prev.W;
        sR = prev.r;
        fprintf('\n  [seed] from theta=%.4g: W=%.4f r=%.4f\n', prev.theta, sW, sR);
    end

    converged = false;
    nev       = nev_seed;

    if USE_BROYDEN
        [Wc, rc, Rnewc, converged, ~, nev_b] = clear_frictionless_broyden( ...
            theta, sW, sR, SEED_J, BIS_TOL, BRO_MAX_EVALS, BRO_BOX);
        nev = nev + nev_b;
    end

    if ~converged
        if USE_BROYDEN
            fprintf('\n  [fallback] Broyden did not converge -- reverting to nested bisection.\n');
        end
        if ~isempty(prev)
            wlo = 0.85*prev.W;
            whi = 1.18*prev.W;
            rlo = max(0.005, prev.r - 0.02);
            rhi = min(0.99*r_ceiling, prev.r + 0.02);
        else
            [wlo, whi, rlo, rhi, probe] = probe_price_bracket(theta, PROBE_W, PROBE_R);
        end
        [Wc, rc, Rnewc, converged, ~, nev2] = clear_frictionless_bisect( ...
            theta, wlo, whi, rlo, rhi, BIS_TOL, BIS_MAX_EVALS);
        nev = nev + nev2;
    end

    x_full = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    evaluate_economy(x_full);
    m  = LAST_MOMENTS;
    g  = gaps_calibration_last;
    AG = aggregates();

    R = AG;
    R.converged = converged;
    R.n_evals   = nev;
    R.gap_Knew  = 100*g(1);
    R.gap_Labor = 100*g(2);
    R.gap_KnKu  = 100*g(3);
    R.LP        = AG.Y / max(AG.L, 1e-12);
    R.top10 = m.top10;
    R.exit = m.exit;
    R.DtoY_m = m.DtoY;
    R.probe = probe;

    if isempty(fr), fr = R; else, fr(end+1) = R; end
    if converged, prev = R; end

    if ~converged
        fprintf('  [WARN] clearing did NOT converge at theta=%.4g -- this theta is unusable.\n', theta);
    end
    if abs(rc - 0.99*r_ceiling) < 1e-4
        fprintf('  [WARN] r pinned at the 1/beta-1 wall -- the capital market did NOT clear.\n');
    end
    if AG.ent < 0.02 || AG.ent > 0.95
        fprintf('  [WARN] ent share = %.3f -- this looks like the COLLAPSED occupational branch\n', AG.ent);
        fprintf('         (the known W knife-edge). Check the probe sign table before believing it.\n');
    end

    fprintf('\n  theta=%-6.4g  W=%.4f r=%.4f Rnew=%.4f\n', theta, Wc, rc, Rnewc);
    fprintf('  Y=%.5f  Y/L=%.5f  TFP=%.5f  Kagg=%.4f  ent=%.4f\n', ...
        AG.Y, R.LP, AG.TFP, AG.Kagg, AG.ent);
    fprintf('  constrained ents = %5.1f%%   mean mu = %.4f   mean k/k* = %.4f\n', ...
        100*AG.constr, AG.mu_mean, AG.kk_mean);
    fprintf('  gaps: Knew %+.1f%%  Labor %+.1f%%   converged=%d  (%d evals)\n', ...
        R.gap_Knew, R.gap_Labor, converged, nev);

    save(CKPT, 'fr', 'FRICT_THETAS', 'PSTAMP');
end

fprintf('\n\n================================================================\n');
fprintf('  NO-USED CF: calibrated vs frictionless\n');
fprintf('================================================================\n\n');
fprintf('  %-10s %8s %8s %9s %9s %8s %8s %8s %5s\n', ...
    'theta','W','r','Y','Y/L','TFP','ent','constr%','conv');
fprintf('  %-10.4g %8.4f %8.4f %9.5f %9.5f %8.5f %8.4f %8.1f %5d   <- CALIBRATED\n', ...
    base.theta, base.W, base.r, base.Y, base.Y/max(base.L,1e-12), base.TFP, ...
    base.ent, 100*base.constr, base.converged);
for j = 1:numel(fr)
    fprintf('  %-10.4g %8.4f %8.4f %9.5f %9.5f %8.5f %8.4f %8.1f %5d\n', ...
        fr(j).theta, fr(j).W, fr(j).r, fr(j).Y, fr(j).LP, fr(j).TFP, ...
        fr(j).ent, 100*fr(j).constr, fr(j).converged);
end

cvg = [fr.converged] > 0;
if ~any(cvg)
    error('nofriction: no frictionless theta converged. Nothing to report.');
end
th_c       = [fr.theta];
th_c(~cvg) = -inf;
[~, ifr] = max(th_c);
fin = fr(ifr);

fprintf('\n----------------------------------------------------------------\n');
fprintf('  HAS THE FRICTIONLESS LIMIT ACTUALLY BEEN REACHED?\n');
fprintf('----------------------------------------------------------------\n');
fprintf('  %-10s %10s %12s %12s %10s\n', 'theta', 'Y', 'dY vs prev', 'constr%', 'mean k/k*');
for j = 1:numel(fr)
    if j == 1
        dYs = '     --';
    else
        dY_j = abs(fr(j).Y - fr(j-1).Y) / max(abs(fr(j-1).Y), 1e-12);
        dYs  = sprintf('%+.3f%%', 100*dY_j);
    end
    fprintf('  %-10.4g %10.5f %12s %11.2f%% %10.4f\n', ...
        fr(j).theta, fr(j).Y, dYs, 100*fr(j).constr, fr(j).kk_mean);
end

ok_limit = false;
if numel(fr) >= 2
    dY_top   = abs(fr(end).Y - fr(end-1).Y) / max(abs(fr(end-1).Y), 1e-12);
    ok_limit = (dY_top < 0.001) && (fin.constr < 0.01);
    fprintf('\n  top two thetas: |dY| = %.3f%%,  constrained at the top = %.2f%%\n', ...
        100*dY_top, 100*fin.constr);
    if ok_limit
        fprintf('  --> LIMIT REACHED. Quote theta = %.4g as the frictionless economy.\n', fin.theta);
    else
        fprintf('  --> NOT YET. Append a larger theta to FRICT_THETAS (only the new one runs).\n');
    end
end

ok_100 = false;
i100   = find(abs([fr.theta] - 100) < 1e-9, 1);
if ~isempty(i100) && i100 < numel(fr)
    dY100  = abs(fr(i100+1).Y - fr(i100).Y) / max(abs(fr(i100).Y), 1e-12);
    ok_100 = (dY100 < 0.001) && (fr(i100).constr < 0.01);
end

gain_Y   = 100*(fin.Y   / base.Y   - 1);
gain_LP  = 100*(fin.LP  / (base.Y/max(base.L,1e-12)) - 1);
gain_TFP = 100*(fin.TFP / base.TFP - 1);

fprintf('\n----------------------------------------------------------------\n');
fprintf('  HEADLINE -- GE GAIN FROM PERFECT CREDIT (CF, NO used capital)\n');
fprintf('  measured at theta = %.4g\n', fin.theta);
fprintf('----------------------------------------------------------------\n');
fprintf('  output Y          %9.5f -> %9.5f    %+7.2f%%\n', base.Y, fin.Y, gain_Y);
fprintf('  labour prod Y/L   %9.5f -> %9.5f    %+7.2f%%   <- Chen''s Table 3 metric\n', ...
    base.Y/max(base.L,1e-12), fin.LP, gain_LP);
fprintf('  TFP               %9.5f -> %9.5f    %+7.2f%%\n', base.TFP, fin.TFP, gain_TFP);
fprintf('  ent share         %9.4f -> %9.4f\n', base.ent, fin.ent);
fprintf('  constrained ents  %8.1f%% -> %8.1f%%\n', 100*base.constr, 100*fin.constr);
fprintf('\n  The thesis says THIS gain must EXCEED the benchmark''s: with no second-hand\n');
fprintf('  market there is more left for a credit market to fix.\n');
fprintf('  Now run, from GE_code/:  >> compare_nofriction\n');
if ~fin.converged
    fprintf('\n  [!] This clear did NOT converge -- do not quote the gain.\n');
end
fprintf('================================================================\n');

nofriction_result = struct( ...
    'economy',      'cf_no_used_capital', ...
    'frict_thetas', FRICT_THETAS, ...
    'rungs',        fr, ...
    'base',         base, ...
    'frictionless', fin, ...
    'gain_Y_pct',   gain_Y, ...
    'gain_LP_pct',  gain_LP, ...
    'gain_TFP_pct', gain_TFP, ...
    'theta100_is_enough', ok_100, ...
    'pass_frictionless',  fin.constr < 0.01, ...
    'params', struct('theta_cal',theta_cal,'eta_p',eta_p,'omega',omega, ...
    'sigma_y',sigma_y,'F',F,'eta_span',eta_span,'beta',beta), ...
    'timestamp', datestr(now));
save('nofriction_result.mat', 'nofriction_result');
fprintf('\n  [saved] nofriction_result.mat\n');
