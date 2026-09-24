clear;
clear global;
clc;
format short;

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
    cf_phi_override R_TARGET WARM_CE WARM_CW

PHI_BASE    = 0.70;
PHI_ROB     = 0.75;
FRICT_THETA = 100000;
CKPT        = 'robustness_phi_ckpt.mat';

script_dir = fileparts(mfilename('fullpath')); cd(script_dir);
code_dir   = fileparts(script_dir);
repo_root  = fileparts(code_dir);
if exist(fullfile(repo_root,'CompEcon','CEtools'),'dir')
    addpath(fullfile(repo_root,'CompEcon','CEtools'),'-end');
elseif exist(fullfile(code_dir,'CompEcon','CEtools'),'dir')
    addpath(fullfile(code_dir,'CompEcon','CEtools'),'-end');
end
if exist('menufun','file')~=2, error('robustness_phi: menufun.m not on path.'); end

alpha                      = 2/3;
delta                      = 0.06;
Lbar                       = 1;
xai                        = 0;
eps_model                  = 1e-6;
sigma_crra                 = 1.5;
gamma                      = 2;
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
nested_clearing_tol        = 0.005;
nested_clearing_tol_Knew   = 0.005;
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

BIS_TOL       = 0.005;
BIS_MAX_EVALS = 160;
BRO_MAX_EVALS = 30;
SEED_J        = [-9.0,-6.1; -3.5,-22.6];
BRO_BOX       = [0.40 0.85 0.020 0.99*r_ceiling];
BIS_W_LO      = 0.44;
BIS_W_HI      = 0.85;
BIS_R_LO      = 0.04;
BIS_R_HI      = 0.99*r_ceiling;

fprintf('\n================================================================\n');
fprintf('  ROBUSTNESS: the new-capital CES weight phi\n');
fprintf('================================================================\n');
fprintf('  benchmark phi = %.2f   robustness phi = %.2f\n', PHI_BASE, PHI_ROB);
fprintf('  calibration HELD FIXED at the benchmark point; prices re-cleared.\n');
fprintf('    theta=%.6f eta_p=%.6f omega=%.6f\n', theta_cal, eta_p, omega);
fprintf('    kappa=%.6f sigma_y=%.6f beta=%.6f\n', kappa, sigma_y, beta);
fprintf('    zeta=%.6f (independent of phi)\n', zeta);
fprintf('  benchmark cleared prices W=%.6f r=%.6f Rnew=%.6f\n', base_W, base_r, base_Rnew);
fprintf('================================================================\n');

VFCACHE = fullfile(script_dir,'vf_warm_cache.mat');
VFBAK   = fullfile(script_dir,'vf_warm_cache_phi070.bak.mat');
if exist(VFCACHE,'file') && ~exist(VFBAK,'file')
    copyfile(VFCACHE, VFBAK);
    fprintf('\n  [vf cache] backed up to %s\n', VFBAK);
end
restore_vf = onCleanup(@() restore_vf_cache(VFBAK, VFCACHE));

PSTAMP = [theta_cal, eta_p, omega, kappa, sigma_y, beta, rho_y, F, PHI_BASE, PHI_ROB];
S = struct();
if exist(CKPT,'file')
    ck = load(CKPT);
    if isfield(ck,'PSTAMP') && numel(ck.PSTAMP)==numel(PSTAMP) ...
            && max(abs(ck.PSTAMP(:)'-PSTAMP)) < 1e-10 && isfield(ck,'S')
        S = ck.S;
        fprintf('\n  [resume] checkpoint loaded: %s\n', strjoin(fieldnames(S)', ', '));
    else
        fprintf('\n  [resume] %s IGNORED -- different parameter point.\n', CKPT);
    end
end

if ~isfield(S,'f070')
    fprintf('\n\n================================================================\n');
    fprintf('  [0] phi = %.2f  FRICTION  (evaluate at the cleared benchmark prices)\n', PHI_BASE);
    fprintf('================================================================\n');
    set_phi(PHI_BASE);
    theta = theta_cal;
    x0 = [base_W; base_r; base_Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    for w = 1:3, evaluate_economy(x0); end
    evaluate_economy(x0);
    S.f070 = pack_point(PHI_BASE, theta, true);
    save(CKPT,'S','PSTAMP'); fprintf('  [ckpt] stage 0 saved.\n');
    report_point('phi=0.70 friction', S.f070);
end

if ~isfield(S,'f075')
    fprintf('\n\n================================================================\n');
    fprintf('  [1] phi = %.2f  FRICTION  (clearing W, r, Rnew)\n', PHI_ROB);
    fprintf('================================================================\n');
    set_phi(PHI_ROB);
    theta = theta_cal;
    W = base_W; r = base_r; Rnew = base_Rnew;
    struct_p = [0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    x0 = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    for w = 1:3, evaluate_economy(x0); end
    [Wc, rc, Rnewc, cvg, ~] = clear_markets(struct_p, W, r, Rnew, ...
        hybrid_nested_price_lb(:), hybrid_nested_price_ub(:));
    xc = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    evaluate_economy(xc);
    S.f075 = pack_point(PHI_ROB, theta, cvg);
    save(CKPT,'S','PSTAMP'); fprintf('  [ckpt] stage 1 saved.\n');
    report_point('phi=0.75 friction', S.f075);
    if ~cvg
        error('robustness_phi: the phi=%.2f friction clear did NOT converge. Do not quote it.', PHI_ROB);
    end
end

if ~isfield(S,'n075')
    fprintf('\n\n================================================================\n');
    fprintf('  [2] phi = %.2f  FRICTIONLESS TWIN  (theta = %.4g)\n', PHI_ROB, FRICT_THETA);
    fprintf('================================================================\n');
    set_phi(PHI_ROB);
    theta = FRICT_THETA;
    sW = S.f075.W;
    sR = S.f075.r;
    fprintf('  [seed] from the phi=%.2f friction point: W=%.4f r=%.4f\n', PHI_ROB, sW, sR);
    [Wc, rc, Rnewc, cvg, ~, nev] = clear_frictionless_broyden( ...
        theta, sW, sR, SEED_J, BIS_TOL, BRO_MAX_EVALS, BRO_BOX);
    if ~cvg
        fprintf('\n  [fallback] Broyden did not converge -- nested bisection.\n');
        wlo = 0.85*sW; whi = 1.30*sW;
        rlo = max(0.005, sR-0.02); rhi = min(0.99*r_ceiling, sR+0.06);
        [Wc, rc, Rnewc, cvg, ~, nev2] = clear_frictionless_bisect( ...
            theta, wlo, whi, rlo, rhi, BIS_TOL, BIS_MAX_EVALS);
        if ~cvg
            fprintf('\n  [fallback 2] widening to the full box.\n');
            [Wc, rc, Rnewc, cvg, ~, nev3] = clear_frictionless_bisect( ...
                theta, BIS_W_LO, BIS_W_HI, BIS_R_LO, BIS_R_HI, BIS_TOL, BIS_MAX_EVALS);
            nev2 = nev2 + nev3;
        end
        nev = nev + nev2;
    end
    xc = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    evaluate_economy(xc);
    S.n075 = pack_point(PHI_ROB, theta, cvg);
    S.n075.n_evals = nev;
    save(CKPT,'S','PSTAMP'); fprintf('  [ckpt] stage 2 saved.\n');
    report_point('phi=0.75 frictionless', S.n075);
    if ~cvg
        error('robustness_phi: the phi=%.2f frictionless clear did NOT converge. Do not quote it.', PHI_ROB);
    end
    if S.n075.constr > 0.01
        fprintf('  [WARN] %.2f%% still constrained at theta=%.4g -- not yet the frictionless limit.\n', ...
            100*S.n075.constr, FRICT_THETA);
    end
end

NFFILE = fullfile(script_dir,'nofriction_result.mat');
if ~exist(NFFILE,'file')
    error(['robustness_phi: %s not found. Run nofriction.m first -- the phi=%.2f ' ...
        'frictionless twin is the comparison this whole check rests on.'], NFFILE, PHI_BASE);
end
NF  = load(NFFILE);
fnm = fieldnames(NF);
NFR = NF.(fnm{1});
if abs(NFR.params.theta_cal - theta_cal) > 1e-6
    error(['robustness_phi: nofriction_result.mat was built at theta=%.6f but the ' ...
        'benchmark is at theta=%.6f. Not the same point.'], NFR.params.theta_cal, theta_cal);
end

b70 = S.f070;  n70 = NFR.frictionless;  base70 = NFR.base;
b75 = S.f075;  n75 = S.n075;

g70 = gains(base70, n70);
g75 = gains(b75,    n75);

fprintf('\n\n');
fprintf('==============================================================================\n');
fprintf('  ROBUSTNESS TO phi, THE NEW-CAPITAL CES WEIGHT\n');
fprintf('  calibration held fixed; prices re-cleared at each phi\n');
fprintf('==============================================================================\n');

fprintf('\n  CLEARED PRICES\n');
fprintf('  %-26s %12s %12s\n', '', sprintf('phi=%.2f',PHI_BASE), sprintf('phi=%.2f',PHI_ROB));
fprintf('  %-26s %12.4f %12.4f\n', 'wage W',                base70.W,    b75.W);
fprintf('  %-26s %12.4f %12.4f\n', 'interest rate r',       base70.r,    b75.r);
fprintf('  %-26s %12.4f %12.4f\n', 'new-capital rent Rnew', base70.Rnew, b75.Rnew);
fprintf('  %-26s %12.4f %12.4f\n', 'used price pu',         base70.pu,   b75.pu);

fprintf('\n  PANEL A: GAIN FROM ELIMINATING THE COLLATERAL CONSTRAINT\n');
fprintf('  (each economy against its own frictionless twin)\n');
fprintf('  %-26s %12s %12s %12s\n', '', sprintf('phi=%.2f',PHI_BASE), sprintf('phi=%.2f',PHI_ROB), 'difference');
fprintf('  %-26s %12.2f %12.2f %+12.2f\n', 'aggregate output (%)',     g70.Y,   g75.Y,   g75.Y-g70.Y);
fprintf('  %-26s %12.2f %12.2f %+12.2f\n', 'output per worker (%)',    g70.LP,  g75.LP,  g75.LP-g70.LP);
fprintf('  %-26s %12.2f %12.2f %+12.2f\n', 'TFP (%)',                  g70.TFP, g75.TFP, g75.TFP-g70.TFP);
fprintf('  %-26s %12.2f %12.2f %+12.2f\n', 'misallocation L (log pt)', g70.L,   g75.L,   g75.L-g70.L);

fprintf('\n  PANEL B: ALLOCATION IN THE CONSTRAINED EQUILIBRIUM\n');
fprintf('  %-26s %12s %12s %12s\n', '', sprintf('phi=%.2f',PHI_BASE), sprintf('phi=%.2f',PHI_ROB), 'difference');
fprintf('  %-26s %12.3f %12.3f %+12.3f\n', 'mean scale K/K*',    base70.kk_mean, b75.kk_mean, b75.kk_mean-base70.kk_mean);
fprintf('  %-26s %12.3f %12.3f %+12.3f\n', 'share constrained',  base70.constr,  b75.constr,  b75.constr-base70.constr);
fprintf('  %-26s %12.3f %12.3f %+12.3f\n', 'aggregate capital',  base70.Kagg,    b75.Kagg,    b75.Kagg-base70.Kagg);
fprintf('  %-26s %12.3f %12.3f %+12.3f\n', 'entrepreneur share', base70.ent,     b75.ent,     b75.ent-base70.ent);

fprintf('\n  TARGETED MOMENTS (nothing is refit at phi=%.2f, so drift is the point)\n', PHI_ROB);
fprintf('  %-26s %10s %12s %12s\n', '', 'data', sprintf('phi=%.2f',PHI_BASE), sprintf('phi=%.2f',PHI_ROB));
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'debt-to-output',        DtoY,                   b70.m.DtoY,      b75.m.DtoY);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'top-decile empl share', top10_emp_share_target, b70.m.top10,     b75.m.top10);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'firm exit rate',        firm_exit_target,       b70.m.exit,      b75.m.exit);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'wage Gini',             wage_gini_target,       b70.m.wGini,     b75.m.wGini);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'real interest rate',    R_TARGET,               b70.r,           b75.r);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'used-share gradient',   beta_size_data,         b70.m.beta_size, b75.m.beta_size);
fprintf('  %-26s %10.4f %12.4f %12.4f\n', 'new-to-used ratio',     KnKu_ratio,             b70.m.KnKu,      b75.m.KnKu);

fprintf('\n  VERDICT\n');
dL = abs(g75.L - g70.L);
rL = 100*dL/max(abs(g70.L),1e-9);
dY = abs(g75.Y - g70.Y);
rY = 100*dY/max(abs(g70.Y),1e-9);
fprintf('    misallocation loss L : %.2f -> %.2f log points  (%.1f%% change)\n', g70.L, g75.L, rL);
fprintf('    output gain          : %.2f -> %.2f percent     (%.1f%% change)\n', g70.Y, g75.Y, rY);
if rL < 10 && rY < 10
    fprintf('    Both headline numbers move less than 10%% when phi goes from %.2f to %.2f.\n', PHI_BASE, PHI_ROB);
    fprintf('    The result is not an artefact of the new-capital weight. Report this row.\n');
else
    fprintf('    A headline number moves more than 10%%. Report the sensitivity honestly\n');
    fprintf('    and justify phi = %.2f on its own evidence rather than on invariance.\n', PHI_BASE);
end

robustness_phi_result = struct( ...
    'phi_base', PHI_BASE, 'phi_rob', PHI_ROB, 'frict_theta', FRICT_THETA, ...
    'f070', b70, 'f075', b75, 'n075', n75, 'n070', n70, 'base070', base70, ...
    'gains070', g70, 'gains075', g75, ...
    'params', struct('theta_cal',theta_cal,'eta_p',eta_p,'omega',omega, ...
        'kappa',kappa,'sigma_y',sigma_y,'beta',beta,'zeta',zeta,'rho_y',rho_y,'F',F));
save('robustness_phi_result.mat','robustness_phi_result');
fprintf('\n  [saved] robustness_phi_result.mat\n');
fprintf('==============================================================================\n');

function set_phi(p)
global phi cf_phi_override WARM_CE WARM_CW

phi             = p;
cf_phi_override = p;
WARM_CE = [];
WARM_CW = [];
fprintf('  [phi] set to %.4f (global and cf_phi_override)\n', p);
end

function P = pack_point(phi_val, theta_val, cvg)
global gaps_calibration_last LAST_MOMENTS

AG          = aggregates();
P           = AG;
P.phi       = phi_val;
P.theta     = theta_val;
P.converged = logical(cvg);
P.LP        = AG.Y / max(AG.L, 1e-12);
P.m         = LAST_MOMENTS;
g           = gaps_calibration_last;
P.gap_Knew  = 100*g(1);
P.gap_Labor = 100*g(2);
P.gap_KnKu  = 100*g(3);
end

function report_point(nm, P)
fprintf('\n  ---- %s ----\n', nm);
fprintf('    W=%.4f r=%.4f Rnew=%.4f qu=%.4f pu=%.4f\n', P.W, P.r, P.Rnew, P.qu, P.pu);
fprintf('    Y=%.5f Y/L=%.5f TFP=%.5f Kagg=%.4f ent=%.4f\n', P.Y, P.LP, P.TFP, P.Kagg, P.ent);
fprintf('    constrained=%.1f%% mean k/k*=%.4f Kn/Ku=%.3f\n', 100*P.constr, P.kk_mean, P.KnKu);
fprintf('    gaps: Knew %+.1f%% Labor %+.1f%% KnKu %+.1f%% converged=%d\n', ...
    P.gap_Knew, P.gap_Labor, P.gap_KnKu, P.converged);
end

function G = gains(f, nf)
lp_f  = f.Y  / max(f.L,  1e-12);
lp_nf = nf.Y / max(nf.L, 1e-12);
G.Y   = 100*(nf.Y/f.Y - 1);
G.LP  = 100*(lp_nf/lp_f - 1);
G.TFP = 100*(nf.TFP/f.TFP - 1);
G.L   = 100*log(nf.TFP/f.TFP);
end

function restore_vf_cache(bak, live)
if exist(bak,'file')
    try
        copyfile(bak, live);
        fprintf('\n  [vf cache] phi=0.70 warm start restored.\n');
    catch
        fprintf('\n  [vf cache] could not restore %s -- delete it by hand.\n', live);
    end
end
end
