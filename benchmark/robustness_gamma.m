clear; clear global; clc;
format short

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
    INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT eta_span ...
    cf_gamma_override R_TARGET WARM_CE WARM_CW nested_clearing_tol_ratio

GAMMAS   = [1.7 2.5];
GAM_BASE = 2.0;
FRICT_THETA = 100000;
CKPT = 'robustness_gamma_ckpt.mat';

script_dir = fileparts(mfilename('fullpath')); cd(script_dir);
code_dir   = fileparts(script_dir);
repo_root  = fileparts(code_dir);
if exist(fullfile(repo_root,'CompEcon','CEtools'),'dir')
    addpath(fullfile(repo_root,'CompEcon','CEtools'),'-end');
elseif exist(fullfile(code_dir,'CompEcon','CEtools'),'dir')
    addpath(fullfile(code_dir,'CompEcon','CEtools'),'-end');
end
if exist('menufun','file')~=2, error('robustness_gamma: menufun.m not on path.'); end

alpha                      = 2/3;
delta                      = 0.06;
Lbar                       = 1;
xai                        = 0;
eps_model                  = 1e-6;
sigma_crra                 = 1.5;
Ny                         = 7;
eta_span                   = 0.79;
eta                        = eta_span;
lambda_u                   = 1.0;
exit_prob                  = 0.0;
psi_u                      = 0;
sigma_perm                 = 0.0;
productivity_weight        = 1.0;
se                         = 0;
rho                        = 0;
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
nested_clearing_quiet      = false;
nested_clearing_max_iter   = 60;
nested_clearing_tol        = 0.002;
nested_clearing_tol_Knew   = 0.002;
nested_clearing_tol_ratio  = 0.005;
calib_moment_weights       = [30; 30; 30; 30; 30; 30];
WGINI_USE_TRIM             = false;

DtoY                   = 0.6000;
KnKu_ratio             = 2.5699;
beta_size_data         = -0.0269;
wage_gini_target       = 0.372;
firm_exit_target       = 0.09782;
top10_emp_share_target = 0.63210;

parameters_benchmark;
theta_cal = theta;
base_W    = W;
base_r    = r;
base_Rnew = Rnew;

r_ceiling              = 1/beta - 1;
hybrid_nested_price_lb = [0.35; 0.005; 0.04];
hybrid_nested_price_ub = [2.00; 0.99*r_ceiling; 0.26];
zeta                   = (delta + kappa) / KnKu_ratio;

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic);
[ygrid, Py, Ey, piy]          = rouwenhorst_y(Ny, rho_y, sigma_y);

BIS_TOL = 0.005;  BIS_MAX_EVALS = 160;  BRO_MAX_EVALS = 30;
SEED_J  = [-9.0,-6.1; -3.5,-22.6];
BRO_BOX = [0.40 0.85 0.020 0.99*r_ceiling];

fprintf('======================================================== \n');
fprintf('SENSITIVITY TO GAMMA \n');
fprintf('======================================================== \n');
fprintf('  gamma swept over %s; the paper uses %.1f\n', mat2str(GAMMAS), GAM_BASE);
fprintf('  calibration held fixed, prices re-cleared at each gamma\n');
fprintf('    theta=%.6f eta_p=%.6f omega=%.6f\n', theta_cal, eta_p, omega);
fprintf('    kappa=%.6f sigma_y=%.6f beta=%.6f phi=%.2f\n', kappa, sigma_y, beta, phi);
fprintf('======================================================== \n');

PSTAMP = [theta_cal, eta_p, omega, kappa, sigma_y, beta, rho_y, F, GAMMAS(:)'];

CACHE = struct('key', {}, 'pt', {});
if exist(CKPT,'file')
    ck = load(CKPT);
    if isfield(ck,'PSTAMP') && numel(ck.PSTAMP)==numel(PSTAMP) ...
            && max(abs(ck.PSTAMP(:)'-PSTAMP)) < 1e-10 && isfield(ck,'CACHE')
        CACHE = ck.CACHE;
        fprintf('\n  [resume] %d clears already done: %s\n', numel(CACHE), ...
            strjoin({CACHE.key}, ', '));
    end
end

S = struct('g', {});
for ig = 1:numel(GAMMAS)
    g = GAMMAS(ig);
    kf = sprintf('g%.2f_fric', g);
    kn = sprintf('g%.2f_nofr', g);

    fprintf('\n======================================================== \n');
    fprintf('GAMMA = %.2f   FRICTION \n', g);
    fprintf('======================================================== \n');
    hit = find(strcmp({CACHE.key}, kf), 1);
    if ~isempty(hit)
        Fr = CACHE(hit).pt;
        fprintf('  [resume] friction clear already done.\n');
    else
        cf_gamma_override = g;  gamma = g;  WARM_CE = [];  WARM_CW = [];
        theta = theta_cal;  W = base_W;  r = base_r;  Rnew = base_Rnew;

        struct_p = [0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
        x0 = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
        for w = 1:3, evaluate_economy(x0); end
        [Wc, rc, Rnewc, cvg] = clear_markets(struct_p, W, r, Rnew, ...
            hybrid_nested_price_lb(:), hybrid_nested_price_ub(:));
        xc = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
        evaluate_economy(xc);
        Fr = pack_point(g, theta, cvg);
        if ~cvg
            error('robustness_gamma: friction clear did not converge at gamma=%.2f.', g);
        end
        CACHE(end+1) = struct('key', kf, 'pt', Fr);
        save(CKPT, 'CACHE', 'PSTAMP');
    end
    report_point(sprintf('gamma=%.2f friction', g), Fr);

    fprintf('\n======================================================== \n');
    fprintf('GAMMA = %.2f   FRICTIONLESS TWIN \n', g);
    fprintf('======================================================== \n');
    hit = find(strcmp({CACHE.key}, kn), 1);
    if ~isempty(hit)
        Nf = CACHE(hit).pt;
        fprintf('  [resume] frictionless clear already done.\n');
    else
        cf_gamma_override = g;  gamma = g;  WARM_CE = [];  WARM_CW = [];
        theta = FRICT_THETA;
        [Wn, rn, Rnewn, cvgn] = clear_frictionless_broyden(theta, Fr.W, Fr.r, ...
            SEED_J, BIS_TOL, BRO_MAX_EVALS, BRO_BOX);
        if ~cvgn
            fprintf('  [fallback] Broyden failed; nested bisection.\n');
            [Wn, rn, Rnewn, cvgn] = clear_frictionless_bisect(theta, ...
                0.85*Fr.W, 1.30*Fr.W, max(0.005, Fr.r-0.02), ...
                min(0.99*r_ceiling, Fr.r+0.06), BIS_TOL, BIS_MAX_EVALS);
        end
        xn = [Wn; rn; Rnewn; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
        evaluate_economy(xn);
        Nf = pack_point(g, theta, cvgn);
        if ~cvgn
            error('robustness_gamma: frictionless clear did not converge at gamma=%.2f.', g);
        end
        CACHE(end+1) = struct('key', kn, 'pt', Nf);
        save(CKPT, 'CACHE', 'PSTAMP');
    end
    report_point(sprintf('gamma=%.2f frictionless', g), Nf);

    S(ig).g = g;  S(ig).fric = Fr;  S(ig).nofr = Nf;
    S(ig).gain_Y = 100*(Nf.Y/Fr.Y - 1);
    S(ig).gain_L = 100*log(Nf.TFP/Fr.TFP);
    save(CKPT, 'S', 'PSTAMP');
    fprintf('  [ckpt] gamma = %.2f saved.\n', g);
end

fprintf('\n\n');
fprintf('============================================================== \n');
fprintf('SENSITIVITY TO GAMMA \n');
fprintf('============================================================== \n');
fprintf('  %-26s', '');
for i = 1:numel(S), fprintf('%12s', sprintf('gamma=%.1f', S(i).g)); end
fprintf('\n  %s\n', repmat('-', 1, 26 + 12*numel(S)));

rows = { 'wage W',                  @(s) s.fric.W
         'interest rate r',         @(s) s.fric.r
         'new-capital rent Rnew',   @(s) s.fric.Rnew
         'used price pu',           @(s) s.fric.pu
         'used share (%)',          @(s) 100*s.fric.ushare_agg
         'share constrained',       @(s) s.fric.constr
         'mean scale K/K*',         @(s) s.fric.kk_mean
         'output gain (%)',         @(s) s.gain_Y
         'misallocation L (log pt)',@(s) s.gain_L };
for j = 1:size(rows,1)
    fprintf('  %-26s', rows{j,1});
    for i = 1:numel(S), fprintf('%12.3f', rows{j,2}(S(i))); end
    fprintf('\n');
end

ib = find(abs([S.g] - GAM_BASE) < 1e-9, 1);
if ~isempty(ib)
    fprintf('\n  Relative to gamma = %.1f:\n', GAM_BASE);
    for i = 1:numel(S)
        if i == ib, continue; end
        fprintf('    gamma=%.1f : L %+.2f log points, output gain %+.2f points\n', ...
            S(i).g, S(i).gain_L - S(ib).gain_L, S(i).gain_Y - S(ib).gain_Y);
    end
    Ls = [S.gain_L];
    fprintf('\n  L is monotone increasing in gamma? %d\n', all(diff(Ls) > 0));
    fprintf('  The paper claims a more substitutable pair raises the reported\n');
    fprintf('  effect, so gamma = %.1f is the conservative end. That claim holds\n', GAM_BASE);
    fprintf('  here only if L at gamma=%.1f exceeds L at gamma=%.1f.\n', ...
        max([S.g]), GAM_BASE);
end

if numel(S) >= 2
    same = abs(S(1).fric.Y - S(2).fric.Y) < 1e-10 && ...
           abs(S(1).fric.TFP - S(2).fric.TFP) < 1e-10;
    if same
        error(['robustness_gamma: gamma=%.2f and gamma=%.2f produced identical ' ...
            'economies, so the override is not biting and the sweep is void.'], ...
            S(1).g, S(2).g);
    end
    fprintf('\n  override check: gamma=%.2f and gamma=%.2f differ (Y by %.2e),\n', ...
        S(1).g, S(2).g, abs(S(1).fric.Y - S(2).fric.Y));
    fprintf('  so cf_gamma_override is biting and the sweep is meaningful.\n');
end

robustness_gamma_result = struct('gammas', GAMMAS, 'base', GAM_BASE, 'S', S, ...
    'params', struct('theta',theta_cal,'eta_p',eta_p,'omega',omega, ...
        'kappa',kappa,'sigma_y',sigma_y,'beta',beta,'phi',phi,'zeta',zeta));
save('robustness_gamma_result.mat','robustness_gamma_result');
fprintf('\n  [saved] robustness_gamma_result.mat\n');
fprintf('============================================================== \n');

function P = pack_point(gam, theta_val, cvg)
global gaps_calibration_last LAST_MOMENTS

AG          = aggregates();
P           = AG;
P.gamma     = gam;
P.theta     = theta_val;
P.converged = logical(cvg);
P.LP        = AG.Y / max(AG.L, 1e-12);
P.m         = LAST_MOMENTS;
P.ushare_agg = (AG.pu*AG.Ku) / max(AG.Kn + AG.pu*AG.Ku, 1e-12);
g           = gaps_calibration_last;
P.gap_Knew  = 100*g(1);
P.gap_Labor = 100*g(2);
P.gap_KnKu  = 100*g(3);
end

function report_point(nm, P)
fprintf('\n  ---- %s ----\n', nm);
fprintf('    W=%.4f r=%.4f Rnew=%.4f qu=%.4f pu=%.4f\n', P.W, P.r, P.Rnew, P.qu, P.pu);
fprintf('    Y=%.5f TFP=%.5f Kagg=%.4f ent=%.4f used share=%.4f\n', ...
    P.Y, P.TFP, P.Kagg, P.ent, P.ushare_agg);
fprintf('    constrained=%.1f%% mean k/k*=%.4f Kn/Ku=%.3f\n', ...
    100*P.constr, P.kk_mean, P.KnKu);
fprintf('    gaps: Knew %+.1f%% Labor %+.1f%% KnKu %+.1f%% converged=%d\n', ...
    P.gap_Knew, P.gap_Labor, P.gap_KnKu, P.converged);
end
