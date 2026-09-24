function [Wc, rc, Rnewc, converged, gaps_out, n_evals] = clear_frictionless_bisect( ...
    theta_in, W_lo, W_hi, r_lo, r_hi, tol, max_evals, verbose)
global beta eta_span eta_p omega kappa sigma_y exit_prob lambda_u ...
    gaps_calibration_last

if nargin < 8 || isempty(verbose), verbose = true; end

r_ceiling = 1/beta - 1;
if r_hi >= r_ceiling
    warning('clear_frictionless_bisect: r_hi=%.4f >= 1/beta-1=%.4f. Clipping to just below.', ...
        r_hi, r_ceiling);
    r_hi = 0.99 * r_ceiling;
end

n_evals  = 0;
W_ABS_LO = 0.20;
W_ABS_HI = 4.00;
W_seed   = [];

    function [gKn, gL, gKK, Rn] = eval_at(WW, rr)
        Rn = rnew_closed_form(rr);
        x_full = [WW; rr; Rn; 0; theta_in; eta_span; eta_p; omega; kappa; ...
            sigma_y;
            exit_prob;
            lambda_u];
        evaluate_economy(x_full);
        g       = gaps_calibration_last;
        gKn     = g(1);
        gL      = g(2);
        gKK     = g(3);
        n_evals = n_evals + 1;
        if verbose
            fprintf('    [bisect %3d] W=%.4f r=%.4f Rnew=%.4f | Kn=%+8.1f%% L=%+8.1f%% KnKu=%+7.1f%%\n', ...
                n_evals, WW, rr, Rn, 100*gKn, 100*gL, 100*gKK);
        end
    end

    function [Wstar, gKn_at, gL_at, gKK_at, ok] = clear_W(rr)
        if isempty(W_seed)
            wl = W_lo;
            wh = W_hi;
        else
            wl = max(W_ABS_LO, 0.80*W_seed);
            wh = min(W_ABS_HI, 1.25*W_seed);
        end

        [~, gL_l] = eval_at(wl, rr);
        n_exp = 0;
        while gL_l < 0 && n_exp < 5 && n_evals < max_evals
            wl = max(W_ABS_LO, 0.70*wl);
            [~, gL_l] = eval_at(wl, rr);
            n_exp = n_exp + 1;
            if wl <= W_ABS_LO + 1e-9, break; end
        end

        [gKn_h, gL_h] = eval_at(wh, rr);
        n_exp = 0;
        while gL_h > 0 && n_exp < 5 && n_evals < max_evals
            wh = min(W_ABS_HI, 1.40*wh);
            [gKn_h, gL_h] = eval_at(wh, rr);
            n_exp = n_exp + 1;
            if wh >= W_ABS_HI - 1e-9, break; end
        end

        ok = (gL_l > 0) && (gL_h < 0);
        if ~ok
            warning(['clear_frictionless_bisect: the labour bracket does NOT straddle at ' ...
                'r=%.4f (gap %+.1f%% at W=%.3f, %+.1f%% at W=%.3f). Returning the ' ...
                'endpoint -- this is NOT a cleared wage.'], ...
                rr, 100*gL_l, wl, 100*gL_h, wh);
            Wstar  = wh;
            gKn_at = gKn_h;
            gL_at  = gL_h;
            gKK_at = NaN;
            return
        end

        Wstar  = wh;
        gKn_at = gKn_h;
        gL_at  = gL_h;
        gKK_at = NaN;
        while (wh - wl) > 5e-4 && n_evals < max_evals
            Wm = 0.5*(wl + wh);
            [gKn_m, gL_m, gKK_m] = eval_at(Wm, rr);
            Wstar  = Wm;
            gKn_at = gKn_m;
            gL_at  = gL_m;
            gKK_at = gKK_m;
            if abs(gL_m) < tol, break; end
            if gL_m > 0, wl = Wm; else, wh = Wm; end
        end
        W_seed = Wstar;
    end

if verbose
    fprintf('\n  [clear_bisect] nested bisection at theta = %.4g   (Chen-style, no damping)\n', theta_in);
    fprintf('    W in [%.3f %.3f],  r in [%.4f %.4f],  tol = %.2f%%,  budget = %d evals\n', ...
        W_lo, W_hi, r_lo, r_hi, 100*tol, max_evals);
    fprintf('    Rnew pinned in closed form; inner W bracket warm-starts across r.\n\n');
end

rl = r_lo;
rh = r_hi;

[W_at_rl, gKn_l] = clear_W(rl);
[W_at_rh, gKn_h] = clear_W(rh);

Wc        = W_at_rh;
rc        = rh;
Rnewc     = rnew_closed_form(rc);
converged = false;
gaps_out  = [gKn_h; NaN; NaN];

if gKn_l < 0
    warning(['clear_frictionless_bisect: the capital gap is already NEGATIVE at r_lo=%.4f ' ...
        '(gap %+.1f%%). Lower r_lo -- the root is below the bracket.'], rl, 100*gKn_l);
    Wc    = W_at_rl;
    rc    = rl;
    Rnewc = rnew_closed_form(rc);
elseif gKn_h > 0
    warning(['clear_frictionless_bisect: the capital gap is STILL POSITIVE at r_hi=%.4f, ' ...
        'which is the 1/beta-1 wall (gap %+.1f%%). Capital demand exceeds savings ' ...
        'at EVERY admissible interest rate, so there is no stationary equilibrium ' ...
        'in the box. This is an economic statement, not a solver failure.'], ...
        rh, 100*gKn_h);
else
    while (rh - rl) > 1e-4 && n_evals < max_evals
        rm = 0.5*(rl + rh);
        [Wm, gKn_m, gL_m, gKK_m] = clear_W(rm);
        Wc       = Wm;
        rc       = rm;
        Rnewc    = rnew_closed_form(rc);
        gaps_out = [gKn_m; gL_m; gKK_m];
        if abs(gKn_m) < tol && abs(gL_m) < tol
            converged = true;
            break
        end
        if gKn_m > 0, rl = rm; else, rh = rm; end
    end
end

if verbose
    if converged
        status = 'CONVERGED';
    else
        status = 'NOT CONVERGED';
    end
    fprintf('\n  [clear_bisect] %s after %d evals: W=%.4f r=%.4f Rnew=%.4f\n', ...
        status, n_evals, Wc, rc, Rnewc);
    fprintf('                 gaps: Knew %+.2f%%  Labor %+.2f%%  KnKu %+.2f%%\n', ...
        100*gaps_out(1), 100*gaps_out(2), 100*gaps_out(3));
    if ~converged && n_evals >= max_evals
        fprintf('                 (ran out of budget -- raise max_evals or tighten the brackets)\n');
    end
end
end
