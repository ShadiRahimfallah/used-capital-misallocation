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
    MKT_DEMAND_SUPPLY INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT ...
    NEST_W NEST_R NEST_RN

if nargin < 1, opts = struct(); end
gv = @(f,d) getdef(opts, f, d);

alpha               = 2/3;
delta               = 0.06;
Lbar                = 1;
xai                 = 0;
eps_model           = 1e-6;
sigma_crra          = 1.5;
phi                 = 0.70;
gamma               = 2;
Ny                  = 7;
eta_span            = 0.79;
eta                 = eta_span;
exit_prob           = 0.0;
psi_u               = 0;
sigma_perm          = 0.0;
productivity_weight = 1.0;
se                  = 0;
rho                 = 0;

theta    = gv('theta',    1.97);
eta_p    = gv('eta_p',    5.70);
omega    = gv('omega',    0.855);
kappa    = gv('kappa',    0.100);
sigma_y  = gv('sigma_y',  0.275);
F        = gv('F',        0.0);
rho_y    = gv('rho_y',    0.95);
beta     = gv('beta',     0.865);
lambda_u = gv('lambda_u', 1.0);

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
nested_clearing_max_iter = gv('max_iter', 25);
nested_clearing_tol      = gv('tol',      0.010);
nested_clearing_tol_Knew = gv('tol_Knew', 0.005);
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
if ~isempty(r_fix)
    r_lb = r_fix;
    r_ub = r_fix;
end
hybrid_nested_price_lb = [0.30; r_lb; 0.03];
hybrid_nested_price_ub = [0.95; r_ub; 0.30];

W    = gv('W0',    0.5681);
r    = gv('r0',    0.0050);
Rnew = gv('Rnew0', 0.0779);
if ~isempty(r_fix), r = r_fix; end
NEST_W  = W;
NEST_R  = r;
NEST_RN = Rnew;

tag = gv('tag', 'eval');
t0  = tic;

fprintf('\n=========================================================\n');
fprintf('  NIGHT_EVAL  tag=%s\n', tag);
fprintf('  th=%.4f etap=%.4f om=%.4f kap=%.5f sgy=%.4f F=%.3f rho_y=%.3f beta=%.4f lam=%.2f\n', ...
    theta, eta_p, omega, kappa, sigma_y, F, rho_y, beta, lambda_u);
if isempty(r_fix)
    fprintf('  closed economy | r box [%.4f, %.4f] | seeds W=%.4f r=%.4f Rn=%.4f | fast=%d\n', ...
        r_lb, r_ub, W, r, Rnew, fast_diagnostic);
else
    fprintf('  SMALL OPEN ECONOMY: r PINNED at %.4f (Knew gap = -NFA/Kn) | seeds W=%.4f Rn=%.4f | fast=%d\n', ...
        r_fix, W, Rnew, fast_diagnostic);
end
fprintf('=========================================================\n');

zeta = (delta + kappa) / KnKu_ratio;
[qu, pu, Rused] = derive_used_prices(r, delta, kappa, Rnew, zeta, psi_u);
[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic);
[ygrid, Py, Ey, piy] = rouwenhorst_y(Ny, rho_y, sigma_y);

struct_p = [0; theta; eta_span; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u];

n_warmup = gv('n_warmup', 2);
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
        save(sprintf('night_%s.mat', tag), 'out');
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

[quc, puc, Rusedc] = derive_used_prices(rc, delta, kappa, Rnewc, zeta, psi_u);
out = struct( ...
    'tag',tag, 'failed',false, ...
    'theta',theta,'eta_p',eta_p,'omega',omega,'kappa',kappa,'sigma_y',sigma_y, ...
    'F',F,'rho_y',rho_y,'beta',beta,'lambda_u',lambda_u,'eta_span',eta_span, ...
    'r_fix', ternary(isempty(r_fix), NaN, r_fix), ...
    'W',Wc,'r',rc,'Rnew',Rnewc,'qu',quc,'pu',puc,'Rused',Rusedc,'zeta',zeta, ...
    'converged',converged,'obj',yobj, ...
    'gap_Knew',100*g(1),'gap_Labor',100*g(2),'gap_KnKu',100*g(3), ...
    'DtoY',m.DtoY,'wGini',m.wGini,'exit',m.exit,'top10',m.top10, ...
    'ent',m.ent,'empl',m.empl,'KnKu',m.KnKu, ...
    'Y',getf(m,'Y'),'Kn',getf(m,'Kn'),'Ku',getf(m,'Ku'), ...
    'KntoY',getf(m,'KntoY'),'KutoY',getf(m,'KutoY'),'KtoY',getf(m,'KtoY'), ...
    'constr',getf(m,'constr'),'constr_top10',getf(m,'constr_top10'), ...
    'wGini_trim',getf(m,'wGini_trim'), 'beta_size',getf(m,'beta_size'), ...
    'diag',dg, 'minutes',toc(t0)/60, 'stamp',datestr(now));

if isfinite(out.Kn) && isfinite(out.Ku) && isfinite(out.Y) && out.Y > 0
    out.KtoY_value = (out.Kn + puc*out.Ku) / out.Y;
else
    out.KtoY_value = NaN;
end

out.NFA_over_Kn = -out.gap_Knew/100;

print_night(out);
save(sprintf('night_%s.mat', tag), 'out');
fprintf('  [saved] night_%s.mat  (%.1f min)\n', tag, out.minutes);
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

function print_night(o)
fprintf('\n---------------- solve_equilibrium result: %s ----------------\n', o.tag);
fprintf('  obj=%.5f  converged=%d  (%.1f min)\n', o.obj, o.converged, o.minutes);
fprintf('  PRICES   W=%.4f  r=%.4f  Rnew=%.4f  Rused=%.4f  qu=%.4f  pu=%.4f\n', ...
    o.W, o.r, o.Rnew, o.Rused, o.qu, o.pu);
fprintf('  GAPS     Knew=%+.2f%%  Labor=%+.2f%%  KnKu=%+.2f%%   [SOE: NFA/Kn=%+.3f]\n', ...
    o.gap_Knew, o.gap_Labor, o.gap_KnKu, o.NFA_over_Kn);
fprintf('  MOMENTS  DtoY=%.4f  wGini=%.4f  exit=%.4f  top10=%.4f  KnKu=%.3f\n', ...
    o.DtoY, o.wGini, o.exit, o.top10, o.KnKu);
fprintf('  AGGS     Y=%.4f  K/Y(phys)=%.3f  K/Y(value)=%.3f  ent=%.4f  empl=%.4f  constr=%.3f\n', ...
    o.Y, o.KtoY, o.KtoY_value, o.ent, o.empl, o.constr);
if isstruct(o.diag) && isfield(o.diag,'ent_top_abil_decile')
    d = o.diag;
    fprintf('  ENTRY    top-ability-decile entrepreneur share = %.4f\n', d.ent_top_abil_decile);
    fprintf('           entrepreneurs in bottom wealth half   = %.4f\n', d.ent_share_bottom_half_wealth);
    fprintf('           poor(bottom50 a) AND top-decile abil, entrepreneur rate = %.4f\n', d.poor_talented_ent_rate);
    fprintf('           median entry threshold a*(top abil decile) = %.4f  (median wealth %.4f)\n', ...
        d.astar_top_abil, d.median_wealth);
end
fprintf('---------------------------------------------------------\n\n');
end
