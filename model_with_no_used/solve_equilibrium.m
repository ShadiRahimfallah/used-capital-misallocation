function out = solve_equilibrium(opts)
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
    INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT NEST_W NEST_R ...
    NEST_RN cf_cf1_zero_used_prices

if nargin < 1, opts = struct(); end
gv = @(f,d) getdef(opts, f, d);

alpha               = 2/3;
delta               = 0.06;
Lbar                = 1;
xai                 = 0;
eps_model           = 1e-6;
sigma_crra          = 1.5;
gamma               = 2;
Ny                  = 7;
eta_span            = 0.79;
eta                 = eta_span;
exit_prob           = 0.0;
sigma_perm          = 0.0;
productivity_weight = 1.0;
se                  = 0;
rho                 = 0;

phi                     = 1.0;
kappa                   = 0.0;
zeta                    = 0.0;
lambda_u                = 0.0;
psi_u                   = 0.0;
qu                      = 0.0;
pu                      = 0.0;
cf_cf1_zero_used_prices = true;

theta   = gv('theta',   1.97);
eta_p   = gv('eta_p',   5.70);
omega   = gv('omega',   0.855);
sigma_y = gv('sigma_y', 0.275);
F       = gv('F',       0.0);
rho_y   = gv('rho_y',   0.95);
beta    = gv('beta',    0.865);

fast_diagnostic            = gv('fast', true);
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
WGINI_USE_TRIM             = false;

nested_clearing_quiet    = false;
nested_clearing_max_iter = gv('max_iter', 30);
nested_clearing_tol      = gv('tol',      0.010);
nested_clearing_tol_Knew = gv('tol_Knew', 0.015);
calib_moment_weights     = [50; 20; 35; 30; 30];

DtoY                   = 0.6000;
KnKu_ratio             = 2.5699;
beta_size_data         = -0.0269;
wage_gini_target       = 0.372;
firm_exit_target       = 0.09782;
top10_emp_share_target = 0.63210;

r_fix = gv('r_fix', []);
r_lb  = gv('r_lb', 0.002);
r_ub  = gv('r_ub', 0.20);
if ~isempty(r_fix), r_lb = r_fix;  r_ub = r_fix;  end
hybrid_nested_price_lb = [0.30; r_lb; 0.03];
hybrid_nested_price_ub = [0.95; r_ub; 0.30];

W = gv('W0',    0.5450);
r = gv('r0',    0.0150);
if ~isempty(r_fix), r = r_fix; end
Rnew    = gv('Rnew0', r + delta);
NEST_W  = W;
NEST_R  = r;
NEST_RN = Rnew;

tag = gv('tag', 'cf');
t0  = tic;

fprintf('\n=========================================================\n');
fprintf('  NIGHT_EVAL_CF (NO USED CAPITAL)  tag=%s\n', tag);
fprintf('  th=%.4f etap=%.4f om=%.4f sgy=%.4f F=%.3f rho_y=%.3f beta=%.4f\n', ...
    theta, eta_p, omega, sigma_y, F, rho_y, beta);
fprintf('  phi=1 kappa=0 zeta=0 lambda_u=0 qu=pu=0 | seeds W=%.4f r=%.4f Rn=%.4f | fast=%d\n', ...
    W, r, Rnew, fast_diagnostic);
if ~isempty(r_fix)
    fprintf('  SMALL OPEN ECONOMY: r PINNED at %.4f\n', r_fix);
end
fprintf('=========================================================\n');

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic);
[ygrid, Py, Ey, piy] = rouwenhorst_y(Ny, rho_y, sigma_y);

struct_p = [0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];

n_warmup = gv('n_warmup', 3);
for wpass = 1:n_warmup
    x_full = [W; r; Rnew; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
    try
        evaluate_economy(x_full);
        g = gaps_calibration_last;
        fprintf('  [warmup %d] gap_L=%+.1f%% gap_Kn=%+.1f%% DtoY=%.3f ent=%.3f\n', ...
            wpass, 100*g(2), 100*g(1), LAST_MOMENTS.DtoY, LAST_MOMENTS.ent);
    catch ME
        fprintf('  [warmup %d FAILED] %s\n', wpass, ME.message);
    end
end

converged = NaN;
Wc        = W;
rc        = r;
Rnewc     = Rnew;
if gv('do_clear', true)
    try
        [Wc, rc, Rnewc, converged] = clear_markets( ...
            struct_p, W, r, Rnew, hybrid_nested_price_lb(:), hybrid_nested_price_ub(:));
    catch ME
        fprintf('  [CLEAR FAILED] %s\n', ME.message);
        out = struct('tag',tag,'failed',true,'err',ME.message);
        save(sprintf('nightcf_%s.mat', tag), 'out');
        return
    end
end
if ~isempty(r_fix), rc = r_fix; end

x_full = [Wc; rc; Rnewc; 0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];
yobj   = evaluate_economy(x_full);
m      = LAST_MOMENTS;
g      = gaps_calibration_last;

dg = struct();
try
    dg = entry_margin_diagnostics();
catch ME
    fprintf('  [entry diag FAILED] %s\n', ME.message);
    dg.err = ME.message;
end

out = struct('tag',tag,'failed',false,'model','cf_noused', ...
    'theta',theta,'eta_p',eta_p,'omega',omega,'kappa',0,'sigma_y',sigma_y, ...
    'F',F,'rho_y',rho_y,'beta',beta,'lambda_u',0,'eta_span',eta_span, ...
    'r_fix', ternary(isempty(r_fix), NaN, r_fix), ...
    'W',Wc,'r',rc,'Rnew',Rnewc,'qu',0,'pu',0,'Rused',NaN,'zeta',0, ...
    'converged',converged,'obj',yobj, ...
    'gap_Knew',100*g(1),'gap_Labor',100*g(2),'gap_KnKu',100*g(3), ...
    'DtoY',m.DtoY,'wGini',m.wGini,'exit',m.exit,'top10',m.top10, ...
    'ent',m.ent,'empl',m.empl,'KnKu',NaN, ...
    'Y',getf(m,'Y'),'Kn',getf(m,'Kn'),'Ku',0, ...
    'KntoY',getf(m,'KntoY'),'KutoY',0,'KtoY',getf(m,'KtoY'), ...
    'constr',getf(m,'constr'),'constr_top10',getf(m,'constr_top10'), ...
    'wGini_trim',getf(m,'wGini_trim'), ...
    'diag',dg,'minutes',toc(t0)/60,'stamp',datestr(now));
out.KtoY_value  = out.KtoY;
out.NFA_over_Kn = -out.gap_Knew/100;

fprintf('\n--- CF result %s: obj=%.5f conv=%d W=%.4f r=%.4f Rnew=%.4f | DtoY=%.4f ent=%.4f Y=%.5f\n', ...
    tag, out.obj, out.converged, out.W, out.r, out.Rnew, out.DtoY, out.ent, out.Y);
save(sprintf('nightcf_%s.mat', tag), 'out');
fprintf('  [saved] nightcf_%s.mat  (%.1f min)\n', tag, out.minutes);
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = getf(s, f)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
