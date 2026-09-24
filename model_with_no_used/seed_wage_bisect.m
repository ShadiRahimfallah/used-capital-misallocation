function [W_seed, gL, gKn, n_evals] = seed_wage_bisect(theta_in, r_fix, W_lo, W_hi, tol, max_evals, verbose)
global eta_span eta_p omega kappa sigma_y exit_prob lambda_u ...
    gaps_calibration_last LAST_MOMENTS

if nargin < 7 || isempty(verbose),   verbose   = true;  end
if nargin < 6 || isempty(max_evals), max_evals = 7;     end
if nargin < 5 || isempty(tol),       tol       = 0.05;  end

n_evals = 0;

    function [gL_, gKn_, ent_] = evalW(WW)
        Rn = rnew_closed_form(r_fix);
        x_full = [WW; r_fix; Rn; 0; theta_in; eta_span; eta_p; omega; kappa; ...
            sigma_y;
            exit_prob;
            lambda_u];
        evaluate_economy(x_full);
        g       = gaps_calibration_last;
        gKn_    = g(1);
        gL_     = g(2);
        ent_    = LAST_MOMENTS.ent;
        n_evals = n_evals + 1;
        if verbose
            fprintf('    [seed %d] W=%.4f r=%.4f | L=%+8.1f%% Kn=%+8.1f%% ent=%.3f\n', ...
                n_evals, WW, r_fix, 100*gL_, 100*gKn_, ent_);
        end
    end

if verbose
    fprintf('\n  [seed_W] coarse 1-D wage bisection at theta=%.4g, r fixed at %.4f\n', theta_in, r_fix);
    fprintf('           W in [%.3f %.3f], tol %.0f%% on the labour gap, budget %d evals.\n', ...
        W_lo, W_hi, 100*tol, max_evals);
    fprintf('           This is a SEED for Broyden, not a solution -- coarse is correct.\n\n');
end

wl = W_lo;
wh = W_hi;
[gL_l, ~] = evalW(wl);
[gL_h, gKn_h, ~] = evalW(wh);

W_seed = wh;
gL     = gL_h;
gKn    = gKn_h;

if gL_l < 0
    warning(['seed_wage_bisect: the labour gap is ALREADY NEGATIVE (%+.1f%%) at the ' ...
        'LOWEST wage %.3f. The clearing wage is BELOW the bracket -- lower W_lo.'], ...
        100*gL_l, wl);
    W_seed = wl;
    gL     = gL_l;
    return
end
if gL_h > 0
    warning(['seed_wage_bisect: the labour gap is STILL POSITIVE (%+.1f%%) at the ' ...
        'HIGHEST wage %.3f. Raise W_hi.'], 100*gL_h, wh);
    return
end

while (wh - wl) > 0.01 && n_evals < max_evals
    Wm = 0.5*(wl + wh);
    [gL_m, gKn_m] = evalW(Wm);
    W_seed = Wm;
    gL     = gL_m;
    gKn    = gKn_m;
    if abs(gL_m) < tol, break; end
    if gL_m > 0, wl = Wm; else, wh = Wm; end
end

if verbose
    fprintf('\n  [seed_W] W_seed = %.4f  (labour gap %+.1f%%, capital gap %+.1f%%)\n', ...
        W_seed, 100*gL, 100*gKn);
    fprintf('           handing (W=%.4f, r=%.4f) to Broyden, which will move BOTH.\n', W_seed, r_fix);
end
end
