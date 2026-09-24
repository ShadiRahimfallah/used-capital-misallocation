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
    INV_NA INV_USE_EIGS INV_CURV MOMENTS_FROM_INVARIANT

FRICT_THETAS = [100000];
CKPT         = 'nofriction_frict_checkpoint.mat';
BASEFILE     = 'nofriction_baseline.mat';

script_dir = fileparts(mfilename('fullpath')); cd(script_dir);
code_dir   = fileparts(script_dir);
repo_root  = fileparts(code_dir);
if exist(fullfile(repo_root,'CompEcon','CEtools'),'dir')
    addpath(fullfile(repo_root,'CompEcon','CEtools'),'-end');
elseif exist(fullfile(code_dir,'CompEcon','CEtools'),'dir')
    addpath(fullfile(code_dir,'CompEcon','CEtools'),'-end');
end
if exist('menufun','file')~=2, error('nofriction: menufun.m not on path.'); end

beta                       = 0.857687323825;
alpha                      = 2/3;
delta                      = 0.06;
Lbar                       = 1;
xai                        = 0;
eps_model                  = 1e-6;
sigma_crra                 = 1.5;
phi                        = 0.70;
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
calib_moment_weights       = [50;
    20;
    35;
    30;
    30];
WGINI_USE_TRIM         = false;
hybrid_nested_price_lb = [0.35;
    0.005;
    0.04];
hybrid_nested_price_ub=[2.00;
    0.99*(1/beta-1);
    0.26];
r_ceiling = 1/beta - 1;

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
r_ceiling                 = 1/beta - 1;
hybrid_nested_price_ub(2) = 0.99*r_ceiling;
zeta = (delta + kappa) / KnKu_ratio;

fprintf('  [nofriction] calibrated point read from parameters_benchmark.m\n');
fprintf('    theta=%.6f eta_p=%.6f omega=%.6f kappa=%.6f sigma_y=%.6f beta=%.6f\n', ...
    theta_cal, eta_p, omega, kappa, sigma_y, beta);
fprintf('    base prices W=%.6f r=%.6f Rnew=%.6f\n', base_W, base_r, base_Rnew);

BIS_TOL       = 0.005;
BIS_MAX_EVALS = 160;
BRO_MAX_EVALS = 30;
USE_BROYDEN   = true; USE_PROBE=true;
SEED_W        = 0.53; SEED_R=0.10;
SEED_J        = [-9.0,-6.1; -3.5,-22.6];
BRO_BOX       = [0.40 0.85 0.020 0.99*r_ceiling];
PROBE_W       = [0.48 0.56 0.64 0.72];
PROBE_R       = [0.05 0.09 0.13];
BIS_W_LO      = 0.44;
BIS_W_HI      = 0.85;
BIS_R_LO      = 0.04;
BIS_R_HI      = 0.99*r_ceiling;

ap_safe = max(eta_p,1.10);
om_safe = max(min(omega,0.999),0.0);
k       = 40;

[egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, false);
[ygrid,Py,Ey,piy]=rouwenhorst_y(Ny,rho_y,sigma_y);

fprintf('================================================================\n');
fprintf('  BENCHMARK (with used capital) -- ELIMINATE THE COLLATERAL CONSTRAINT\n');
fprintf('================================================================\n');
fprintf('  final point: theta_cal=%.4f eta_p=%.4f omega=%.4f kappa=%.4f sigma_y=%.4f F=%.4f rho_y=%.3f\n', ...
    theta_cal, eta_p, omega, kappa, sigma_y, F, rho_y);
fprintf('  NOT recalibrated; only prices re-cleared. frictionless thetas: %s\n', mat2str(FRICT_THETAS));
fprintf('  r wall 1/beta-1 = %.4f (no stationary eq above it)\n', r_ceiling);
fprintf('================================================================\n');

PSTAMP = [theta_cal, eta_p, omega, kappa, sigma_y, beta, rho_y, F];

base = [];
if exist(BASEFILE,'file')
    bb = load(BASEFILE);
    if isfield(bb,'PSTAMP') && numel(bb.PSTAMP)==numel(PSTAMP) && max(abs(bb.PSTAMP(:)'-PSTAMP)) < 1e-10
        base = bb.base;
        fprintf('\n  [baseline] loaded from %s (parameter stamp matches)\n', BASEFILE);
    else
        fprintf('\n  [baseline] %s IGNORED -- it was built at different parameters. Re-solving.\n', BASEFILE);
    end
end
if isempty(base)
    fprintf('\n  [baseline] clearing theta=%.4f once (seeded at its cleared prices) and caching.\n', theta_cal);
    theta    = theta_cal;
    W        = base_W;
    r        = base_r;
    Rnew     = base_Rnew;
    struct_p = [0;
        theta;
        eta_span;
        eta_p;
        omega;
        kappa;
        sigma_y;
        exit_prob;
        lambda_u];
    x_full=[W;
        r;
        Rnew;
        0;
        theta;
        eta_span;
        eta_p;
        omega;
        kappa;
        sigma_y;
        exit_prob;
        lambda_u];
    for w=1:3, evaluate_economy(x_full); end
    [Wc,rc,Rnewc,cvg,~]=clear_markets(struct_p,W,r,Rnew,hybrid_nested_price_lb(:),hybrid_nested_price_ub(:));
    x_full=[Wc;
        rc;
        Rnewc;
        0;
        theta;
        eta_span;
        eta_p;
        omega;
        kappa;
        sigma_y;
        exit_prob;
        lambda_u];
    evaluate_economy(x_full);
    gb   = gaps_calibration_last;
    base = aggregates();
    base.converged=cvg;
    base.gap_Knew=100*gb(1);
    base.gap_Labor=100*gb(2);
    base.gap_KnKu=100*gb(3);
    base.LP=base.Y/max(base.L,1e-12);
    if ~cvg, error('nofriction: the BASELINE clear did not converge -- fix before proceeding.'); end
    save(BASEFILE,'base','PSTAMP'); fprintf('  [baseline] cleared and saved to %s\n', BASEFILE);
end
fprintf('  [baseline] theta=%.4f W=%.4f r=%.4f Rnew=%.4f | Y=%.5f TFP=%.5f ent=%.4f constr=%.1f%% k/k*=%.4f\n', ...
    base.theta, base.W, base.r, base.Rnew, base.Y, base.TFP, base.ent, 100*base.constr, base.kk_mean);

fr   = struct([]);
j0   = 1;
prev = [];
if exist(CKPT,'file')
    ck=load(CKPT);
    if isfield(ck,'PSTAMP') && (numel(ck.PSTAMP)~=numel(PSTAMP) || max(abs(ck.PSTAMP(:)'-PSTAMP)) >= 1e-10)
        fprintf('\n  [resume] %s IGNORED -- built at different parameters. Re-solving.\n', CKPT);
        ck = struct();
    elseif ~isfield(ck,'PSTAMP')
        fprintf('\n  [resume] %s IGNORED -- no parameter stamp (pre-dates the guard). Re-solving.\n', CKPT);
        ck = struct();
    end
    if isfield(ck,'fr') && isfield(ck,'FRICT_THETAS')
        nd=numel(ck.FRICT_THETAS);
        if nd<=numel(FRICT_THETAS) && isequal(ck.FRICT_THETAS(:)',FRICT_THETAS(1:nd))
            fr=ck.fr; j0=numel(fr)+1; ip=find([fr.converged],1,'last');
            if ~isempty(ip), prev=fr(ip); end
            fprintf('\n  [resume] %d frictionless theta(s) done; continuing.\n', numel(fr));
        end
    end
end

for j=j0:numel(FRICT_THETAS)
    theta=FRICT_THETAS(j);
    fprintf('\n================================================================\n');
    fprintf('  FRICTIONLESS CLEAR:  theta = %.4g\n', theta);
    fprintf('================================================================\n');
    probe=[];
    if isempty(prev), sW=SEED_W; sR=SEED_R; else, sW=prev.W; sR=prev.r; fprintf('  [seed] from theta=%.4g: W=%.4f r=%.4f\n', prev.theta, sW, sR); end
    converged = false;
    nev       = 0;
    if USE_BROYDEN
        [Wc,rc,Rnewc,converged,~,nev]=clear_frictionless_broyden(theta,sW,sR,SEED_J,BIS_TOL,BRO_MAX_EVALS,BRO_BOX);
    end
    if ~converged
        if USE_BROYDEN, fprintf('\n  [fallback] Broyden did not converge -- probe + nested bisection.\n'); end
        if USE_PROBE && isempty(prev)
            [wlo,whi,rlo,rhi,probe]=probe_price_bracket(theta,PROBE_W,PROBE_R);
        elseif isempty(prev)
            wlo = BIS_W_LO;
            whi = BIS_W_HI;
            rlo = BIS_R_LO;
            rhi = BIS_R_HI;
        else
            wlo = 0.85*prev.W;
            whi = 1.18*prev.W;
            rlo = max(0.005,prev.r-0.02);
            rhi = min(0.99*r_ceiling,prev.r+0.02);
        end
        [Wc,rc,Rnewc,converged,~,nev2]=clear_frictionless_bisect(theta,wlo,whi,rlo,rhi,BIS_TOL,BIS_MAX_EVALS);
        nev=nev+nev2;
    end

    x_full=[Wc;
        rc;
        Rnewc;
        0;
        theta;
        eta_span;
        eta_p;
        omega;
        kappa;
        sigma_y;
        exit_prob;
        lambda_u];
    evaluate_economy(x_full);
    m  = LAST_MOMENTS;
    g  = gaps_calibration_last;
    AG = aggregates();
    R  = AG;
    R.converged=converged;
    R.n_evals=nev;
    R.gap_Knew=100*g(1);
    R.gap_Labor=100*g(2);
    R.gap_KnKu=100*g(3);
    R.LP=AG.Y/max(AG.L,1e-12);
    R.top10=m.top10;
    R.exit=m.exit;
    R.DtoY_m=m.DtoY;
    R.probe=probe;
    if isempty(fr), fr=R; else, fr(end+1)=R; end
    if converged, prev=R; end
    if ~converged, fprintf('  [WARN] did NOT converge at theta=%.4g -- unusable.\n', theta); end
    if abs(R.gap_KnKu)>2.0
        fprintf('  [WARN] KnKu gap %+.1f%%: closed-form Rnew assumes nobody constrained; theta=%.4g still binds for some.\n', R.gap_KnKu, theta);
    end
    fprintf('\n  theta=%-6.4g W=%.4f r=%.4f Rnew=%.4f qu=%.4f pu=%.4f\n', theta, Wc, rc, Rnewc, AG.qu, AG.pu);
    fprintf('  Y=%.5f Y/L=%.5f TFP=%.5f Kagg=%.4f ent=%.4f\n', AG.Y, R.LP, AG.TFP, AG.Kagg, AG.ent);
    fprintf('  constrained ents=%.1f%% mean mu=%.4f mean k/k*=%.4f\n', 100*AG.constr, AG.mu_mean, AG.kk_mean);
    fprintf('  gaps: Knew %+.1f%% Labor %+.1f%% KnKu %+.1f%% converged=%d (%d evals)\n', R.gap_Knew, R.gap_Labor, R.gap_KnKu, converged, nev);
    save(CKPT,'fr','FRICT_THETAS','PSTAMP');
end

fprintf('\n\n================================================================\n');
fprintf('  BENCHMARK: calibrated vs frictionless\n');
fprintf('================================================================\n\n');
fprintf('  %-10s %8s %8s %9s %9s %8s %8s %8s %5s\n','theta','W','r','Y','Y/L','TFP','ent','constr%','conv');
fprintf('  %-10.4g %8.4f %8.4f %9.5f %9.5f %8.5f %8.4f %8.1f %5d   <- CALIBRATED\n', ...
    base.theta, base.W, base.r, base.Y, base.Y/max(base.L,1e-12), base.TFP, base.ent, 100*base.constr, base.converged);
for j=1:numel(fr)
    fprintf('  %-10.4g %8.4f %8.4f %9.5f %9.5f %8.5f %8.4f %8.1f %5d\n', ...
        fr(j).theta, fr(j).W, fr(j).r, fr(j).Y, fr(j).LP, fr(j).TFP, fr(j).ent, 100*fr(j).constr, fr(j).converged);
end

cvg=[fr.converged]>0;
if ~any(cvg), error('nofriction: no frictionless theta converged.'); end
th_c       = [fr.theta];
th_c(~cvg) = -inf;
[~,ifr]=max(th_c);
fin=fr(ifr);

fprintf('\n----------------------------------------------------------------\n');
fprintf('  FRICTIONLESS LIMIT REACHED?  (need constr~0 AND Y stopped moving)\n');
fprintf('----------------------------------------------------------------\n');
fprintf('  %-10s %10s %12s %12s %10s\n','theta','Y','dY vs prev','constr%','mean k/k*');
for j=1:numel(fr)
    if j==1, dYs='     --'; else, dYs=sprintf('%+.3f%%',100*abs(fr(j).Y-fr(j-1).Y)/max(abs(fr(j-1).Y),1e-12)); end
    fprintf('  %-10.4g %10.5f %12s %11.2f%% %10.4f\n', fr(j).theta, fr(j).Y, dYs, 100*fr(j).constr, fr(j).kk_mean);
end
ok_limit=false;
if numel(fr)>=2
    dY_top   = abs(fr(end).Y-fr(end-1).Y)/max(abs(fr(end-1).Y),1e-12);
    ok_limit = (dY_top<0.001)&&(fin.constr<0.01);
    fprintf('\n  top two thetas: |dY|=%.3f%%, constrained at the top=%.2f%%\n', 100*dY_top, 100*fin.constr);
    if ok_limit, fprintf('  --> LIMIT REACHED. Quote theta=%.4g as the frictionless economy.\n', fin.theta);
    else, fprintf('  --> NOT YET. Append a larger theta (only the new one runs).\n'); end
else
    ok_limit=(fin.constr<0.01);
    if ok_limit, fprintf('\n  single theta=%.4g: constrained=%.2f%% ~ 0 -> frictionless economy. (Add a 2nd larger theta to auto-check Y-invariance.)\n', fin.theta, 100*fin.constr);
    else, fprintf('\n  single theta=%.4g: constrained=%.2f%% still > 0 -> raise theta.\n', fin.theta, 100*fin.constr); end
end

gain_Y   = 100*(fin.Y/base.Y-1);
gain_LP  = 100*(fin.LP/(base.Y/max(base.L,1e-12))-1);
gain_TFP = 100*(fin.TFP/base.TFP-1);
fprintf('\n----------------------------------------------------------------\n');
fprintf('  HEADLINE -- GE GAIN FROM PERFECT CREDIT (benchmark, WITH used capital)\n');
fprintf('  at theta=%.4g\n', fin.theta);
fprintf('----------------------------------------------------------------\n');
fprintf('  output Y         %9.5f -> %9.5f   %+7.2f%%\n', base.Y, fin.Y, gain_Y);
fprintf('  labour prod Y/L  %9.5f -> %9.5f   %+7.2f%%\n', base.Y/max(base.L,1e-12), fin.LP, gain_LP);
fprintf('  TFP              %9.5f -> %9.5f   %+7.2f%%\n', base.TFP, fin.TFP, gain_TFP);
fprintf('  r (deposit rate) %9.4f -> %9.4f   (glut dissipates -> r rises)\n', base.r, fin.r);
fprintf('  ent share        %9.4f -> %9.4f\n', base.ent, fin.ent);
fprintf('  constrained ents %8.1f%% -> %8.1f%%\n', 100*base.constr, 100*fin.constr);
fprintf('  mean k/k*        %9.4f -> %9.4f   (1.0 = unconstrained scale)\n', base.kk_mean, fin.kk_mean);
fprintf('\n  THESIS: compare this to cf_noused/nofriction (perfect credit WITHOUT used capital).\n');
fprintf('  This gain should be the SMALLER of the two -- used capital already does part of a\n');
fprintf('  credit market''s job, so perfect credit has less left to fix.\n');
if ~fin.converged, fprintf('\n  [!] This clear did NOT converge -- do not quote the gain.\n'); end
fprintf('================================================================\n');

nofriction_result=struct('economy','benchmark_used_capital','frict_thetas',FRICT_THETAS, ...
    'rungs',fr,'base',base,'frictionless',fin,'gain_Y_pct',gain_Y,'gain_LP_pct',gain_LP,'gain_TFP_pct',gain_TFP, ...
    'pass_frictionless',fin.constr<0.01, ...
    'params',struct('theta_cal',theta_cal,'eta_p',eta_p,'omega',omega,'kappa',kappa,'sigma_y',sigma_y, ...
    'F',F,'eta_span',eta_span,'lambda_u',lambda_u,'rho_y',rho_y,'beta',beta));
save('nofriction_result.mat','nofriction_result');
fprintf('\n  [saved] nofriction_result.mat\n');
