clear;
clc;
here = fileparts(mfilename('fullpath'));
cd(here);

diary('run_cf_ge.log'); diary on;
fprintf('==== run_cf_ge  (no used capital, recalibrated, r targeted)  %s ====\n\n', datestr(now));

parameters_cf_calib;

opts = struct( ...
    'theta',    theta,   ...
    'eta_p',    eta_p,   ...
    'omega',    omega,   ...
    'sigma_y',  sigma_y, ...
    'beta',     beta,    ...
    'F',        F,       ...
    'rho_y',    rho_y,   ...
    'do_clear', true,    ...
    'fast',     false,   ...
    'tag',      'model_with_no_used');

opts.W0    = W;
opts.r0    = r;
opts.Rnew0 = Rnew;

fprintf('  theta=%.4f eta_p=%.4f omega=%.4f sigma_y=%.4f beta=%.4f (FREE - search kept it here)\n', ...
    opts.theta, opts.eta_p, opts.omega, opts.sigma_y, opts.beta);
fprintf('  used capital removed: kappa=0, lambda_u=0, phi=1, qu=pu=0 (set inside solve_equilibrium)\n\n');

o = solve_equilibrium(opts);

TGT = struct('DtoY',0.6000,'wGini',0.372,'exit',0.09782,'top10',0.6321);
mm  = {'DtoY','wGini','exit','top10'};

fprintf('\n================= CF fully recalibrated, beta free =================\n');
if isfield(o,'failed') && o.failed
    fprintf('  FAILED: %s\n', o.err);
else
    fprintf('  converged = %d\n', o.converged);
    fprintf('  PRICES  W=%.4f r=%.4f Rnew=%.4f\n', o.W, o.r, o.Rnew);
    fprintf('  GAPS    Knew=%+.2f%%  Labor=%+.2f%%\n', o.gap_Knew, o.gap_Labor);
    fprintf('  Y       = %.6f\n', o.Y);
    fprintf('  ent     = %.4f   empl = %.4f   constr = %.4f\n', o.ent, o.empl, o.constr);
    fprintf('\n  %-8s %10s %10s %9s\n', 'moment', 'model', 'target', 'gap%');
    for j = 1:numel(mm)
        mv = o.(mm{j});
        tv = TGT.(mm{j});
        fprintf('  %-8s %10.4f %10.4f %+8.1f%%\n', mm{j}, mv, tv, 100*(mv-tv)/abs(tv));
    end
    fprintf('\n  r (TARGETED in this variant) = %.4f  vs 0.0310\n', o.r);
end

save('run_cf_ge_result.mat', 'o', 'opts');
fprintf('\n==== run_cf_ge DONE %s ====\n', datestr(now));
diary off;
