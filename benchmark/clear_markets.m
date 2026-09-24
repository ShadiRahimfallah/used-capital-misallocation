function [W, r, Rnew, converged, gaps_out] = clear_markets(struct_params, W0, r0, Rnew0, lb_p, ub_p)
global gaps_calibration_last calib_mode delta hybrid_inner_clearing_skip ...
    nested_clearing_quiet nested_clearing_max_iter nested_clearing_tol ...
    nested_clearing_alpha_init KnKu_ratio nested_clearing_tol_ratio

if isempty(hybrid_inner_clearing_skip), hybrid_inner_clearing_skip = false; end

if isempty(KnKu_ratio)
    KnKu_ratio_local = 2.5699;
else
    KnKu_ratio_local = KnKu_ratio;
end

if isempty(delta)
    delta_local = 0.06;
    fprintf('  [clear_markets] WARNING: global delta empty, using fallback %.2f\n', delta_local);
else
    delta_local = delta;
end

szp = numel(struct_params(:));
if szp == 8
    struct_exp = [struct_params(1:2); 0.79; struct_params(3:8)];
elseif szp == 9
    struct_exp = struct_params(:);
else
    error(['clear_markets: struct_params must have 8 (legacy, no eta) or 9 ', ...
        '[reserved_slot; theta; eta; eta_p; omega; kappa; sigma_y; exit_prob; lambda_u]. Got %d.'], szp);
end

quiet = (~isempty(nested_clearing_quiet)) && logical(nested_clearing_quiet(1));

if isempty(nested_clearing_max_iter) || nested_clearing_max_iter < 1
    max_evals = 40;
else
    max_evals = nested_clearing_max_iter;
end

if isempty(nested_clearing_tol) || nested_clearing_tol <= 0
    tol = 0.01;
else
    tol = nested_clearing_tol;
end

if ~isempty(nested_clearing_tol_ratio) && nested_clearing_tol_ratio > 0
    tol_ratio = nested_clearing_tol_ratio;
else
    tol_ratio = max(4*tol, 0.02);
end

if ~isempty(nested_clearing_alpha_init) && nested_clearing_alpha_init > 0
    alpha_init = nested_clearing_alpha_init;
else
    alpha_init = 0.15;
end

alpha_damp_min = 0.03;
stall_tol      = 0.002;
stall_max      = 14;
max_inner      = 5;
eps_s          = -0.42;

gap_clip       = 0.30;
flip_gap       = 0.60;
flip_bail      = 3;
climb_bail     = 12;
climb_step     = 4e-3;
cliff_pad      = 5e-4;
w_climb_step   = 0.04;

kappa_p  = struct_exp(6);
zeta_loc = (delta_local + kappa_p) / KnKu_ratio_local;

W          = W0;
r          = r0;
s          = Rnew0 - r0 - delta_local;
s          = min(max(s, 1e-4), zeta_loc - 1e-4);
alpha_damp = alpha_init;

best_obj    = inf;
best_W      = W0;
best_r      = r0;
best_Rnew   = Rnew0;
best_gaps   = [];
prev_gaps   = [];
prev_mkt    = inf;
stall_count = 0;
converged   = false;
n_evals     = 0;
outer       = 0;
n_flips     = 0;
r_cliff     = -inf;
r_prev      = r0;
r_flip_lo   = -inf;
W_cliff     = -inf;

if ~quiet
    fprintf('\n  [clear_markets] nested clearing of (W, r, Rnew): budget %d evals, tol=%.1f%% (ratio %.1f%%)\n', ...
        max_evals, 100*tol, 100*tol_ratio);
end

while n_evals < max_evals
    outer = outer + 1;

    s_lo       = 1e-4;
    s_hi       = zeta_loc - 1e-4;
    s_prev     = NaN;
    gap_prev   = NaN;
    g          = [];
    gap_KnKu   = NaN;
    degenerate = false;

    for it_s = 1:max_inner
        if n_evals >= max_evals, break; end

        Rnew_try = r + delta_local + s;
        Rnew_try = max(min(Rnew_try, ub_p(3)), lb_p(3));
        s        = Rnew_try - r - delta_local;

        g       = eval_gaps(W, r, Rnew_try, struct_exp);
        n_evals = n_evals + 1;

        if isempty(g) || numel(g) < 3
            if ~quiet, fprintf('    outer %2d: degenerate evaluation (penalty), abandoning point\n', outer); end
            degenerate = true;
            break;
        end

        gap_KnKu = g(3);
        if ~quiet
            [qu_p, pu_p, Rused_p] = derive_used_prices(r, delta_local, kappa_p, Rnew_try, zeta_loc, 0);
            fprintf('      [s %d.%d] s=%.5f pu=%.4f qu=%.4f Rused=%.4f | KnKu=%+.1f%%\n', ...
                outer, it_s, s, pu_p, qu_p, Rused_p, 100*gap_KnKu);
        end
        if abs(gap_KnKu) < tol_ratio, break; end

        if it_s > 1 && isfinite(gap_prev) && abs(log(s/s_prev)) > 1e-8
            eps_est = log((1+gap_KnKu)/(1+gap_prev)) / log(s/s_prev);
            if isfinite(eps_est) && eps_est < -0.05, eps_s = eps_est; end
        end
        s_prev   = s;
        gap_prev = gap_KnKu;

        if gap_KnKu > 0
            s_lo = max(s_lo, s);
        else
            s_hi = min(s_hi, s);
        end

        s = exp(log(s) + log(1 + gap_KnKu) / (-eps_s));
        if s <= s_lo || s >= s_hi, s = sqrt(s_lo * s_hi); end
        s = min(max(s, s_lo*1.0001), s_hi*0.9999);
    end
    if isempty(g), break; end

    gap_Knew  = g(1);
    gap_Labor = g(2);
    gap_KnKu  = g(3);
    Rnew      = r + delta_local + s;

    flipped_basin = degenerate || max(abs([gap_Knew, gap_Labor])) > flip_gap;
    if flipped_basin
        n_flips       = n_flips + 1;
        worker_flood  = ~degenerate && gap_Labor < -flip_gap;
        ent_explosion = ~degenerate && gap_Labor >  flip_gap;

        if ~ent_explosion
            r_flip_lo = max(r_flip_lo, r);
        end

        if worker_flood
            r_cliff = max(r_cliff, r);
            if ~quiet
                fprintf('    outer %2d: ENT-COLLAPSE at r=%.4f (L=%+.0f%%) -- r climbed above %.4f for this eval\n', ...
                    outer, r, 100*gap_Labor, r_cliff + climb_step);
            end
        elseif ent_explosion
            W_cliff = max(W_cliff, W);
            if ~quiet
                fprintf('    outer %2d: ENT-EXPLOSION at W=%.4f (L=%+.0f%%) -- W climbed above %.4f for this eval [flip %d]\n', ...
                    outer, W, 100*gap_Labor, min(W_cliff*(1+w_climb_step), ub_p(1)), n_flips);
            end
        elseif ~quiet
            fprintf('    outer %2d: BASIN FLIP (Kn=%+.0f%% L=%+.0f%%) -- rollback + halve damp [flip %d]\n', ...
                outer, 100*gap_Knew, 100*gap_Labor, n_flips);
        end

        climb_mode = isempty(best_gaps) && (worker_flood || ent_explosion);
        this_bail  = flip_bail;
        if climb_mode, this_bail = climb_bail; end

        if n_flips >= this_bail
            if ~quiet
                if climb_mode
                    fprintf('    CLIFF UNREACHABLE: %d climbs, r=%.4f W=%.4f -- returning best (mkt=%.3f) after %d evals\n', ...
                        n_flips, r, W, best_obj, n_evals);
                else
                    fprintf('    BISTABLE point: %d basin flips -- returning best (mkt=%.3f) after %d evals\n', ...
                        n_flips, best_obj, n_evals);
                end
            end
            break;
        end

        if ~isempty(best_gaps)
            W = best_W;
            r = best_r;
            s = best_Rnew - best_r - delta_local;
        else
            W = W0;
            r = r0;
            s = Rnew0 - r0 - delta_local;
        end

        pad = cliff_pad;
        if climb_mode, pad = climb_step; end
        r   = min(max(r, r_cliff + pad), ub_p(2));

        if isfinite(W_cliff)
            W = min(max(W, W_cliff * (1 + w_climb_step)), ub_p(1));
        end

        s          = min(max(s, 1e-4), zeta_loc - 1e-4);
        alpha_damp = max(alpha_damp * 0.5, alpha_damp_min);
        prev_gaps  = [];
        r_prev     = r;
        continue;
    end

    mkt_obj = sqrt(gap_Knew^2 + gap_Labor^2 + 0.25*gap_KnKu^2);
    if mkt_obj < best_obj
        best_obj  = mkt_obj;
        best_W    = W;
        best_r    = r;
        best_Rnew = Rnew;
        best_gaps = g;
    end

    if ~quiet
        [qu_p, pu_p, ~] = derive_used_prices(r, delta_local, kappa_p, Rnew, zeta_loc, 0);
        fprintf(['    iter %2d: W=%.4f r=%.4f Rnew=%.4f s=%.4f qu=%.4f pu=%.4f | ' ...
            'Kn=%+.1f%% L=%+.1f%% KnKu=%+.1f%% | mkt=%.3f (evals %d/%d)\n'], ...
            outer, W, r, Rnew, s, qu_p, pu_p, ...
            100*gap_Knew, 100*gap_Labor, 100*gap_KnKu, mkt_obj, n_evals, max_evals);
    end

    std_conv = abs(gap_Knew) < tol && abs(gap_Labor) < tol && abs(gap_KnKu) < tol_ratio;

    KN_CLIFF_BAND = 0.005;
    cliff_pinned  = (r_cliff > 0) && (gap_Knew < 0) && (abs(gap_Knew) < KN_CLIFF_BAND) && ...
        abs(gap_Labor) < tol && abs(gap_KnKu) < tol_ratio;

    if std_conv || cliff_pinned
        converged = true;
        if ~quiet
            if cliff_pinned && ~std_conv
                fprintf('    CLIFF-PINNED convergence at outer %d: r=%.4f at cliff floor, Labor+KnKu clear, Knew=%+.1f%% structural (%d evals)\n', ...
                    outer, r, 100*gap_Knew, n_evals);
            else
                fprintf('    CONVERGED at outer %d (all gaps in tol, %d evals)\n', outer, n_evals);
            end
        end
        break;
    end

    if outer > 3
        improvement = (prev_mkt - mkt_obj) / max(prev_mkt, 0.01);
        if improvement < stall_tol
            stall_count = stall_count + 1;
        else
            stall_count = 0;
        end
        if stall_count >= stall_max
            if ~quiet
                fprintf('    STALLED at outer %d (improv < %.1f%% for %d iters, best mkt=%.3f)\n', ...
                    outer, 100*stall_tol, stall_max, best_obj);
            end
            break;
        end
    end
    prev_mkt = mkt_obj;

    if ~isempty(prev_gaps)
        sflip = (sign(gap_Knew) ~= sign(prev_gaps(1))) || ...
            (sign(gap_Labor) ~= sign(prev_gaps(2)));
        if sflip
            alpha_damp = max(alpha_damp * 0.65, alpha_damp_min);
        else
            alpha_damp = min(alpha_damp * 1.10, 0.30);
        end
    end

    use_secant = false;
    r_secant   = r;
    if ~isempty(prev_gaps) && abs(gap_Labor) < 1.25*tol && abs(gap_Knew) > 1.5*tol ...
            && abs(r - r_prev) > 1e-7
        slope_Kn = (gap_Knew - prev_gaps(1)) / (r - r_prev);
        if isfinite(slope_Kn) && abs(slope_Kn) > 1e-3
            dr_sec = -0.8 * gap_Knew / slope_Kn;
            dr_sec = max(min(dr_sec, 0.005), -0.005);
            if sign(dr_sec) == sign(gap_Knew) || gap_Knew == 0
                r_secant   = r + dr_sec;
                use_secant = true;
                if isfinite(r_flip_lo) && r_secant <= r_flip_lo + cliff_pad
                    r_secant = r_flip_lo + cliff_pad;
                    if r_secant >= r, use_secant = false; end
                end
            end
        end
    end
    prev_gaps = [gap_Knew; gap_Labor];

    gK     = max(min(gap_Knew,  gap_clip), -gap_clip);
    gL     = max(min(gap_Labor, gap_clip), -gap_clip);
    r_prev = r;
    r_new  = r * (1 + alpha_damp * gK);
    if use_secant
        r_new = r_secant;
        if ~quiet
            fprintf('    [secant] r step %.4f -> %.4f (slope dKn/dr=%.2f, gap=%+.1f%%)\n', ...
                r, r_new, slope_Kn, 100*gap_Knew);
        end
    end
    W_new = W * (1 + alpha_damp * gL);

    W_new = max(min(W_new, ub_p(1)), lb_p(1));
    r_new = max(min(r_new, ub_p(2)), lb_p(2));
    r_new = max(r_new, r_cliff + cliff_pad);

    W = W_new;
    r = r_new;
    s = min(max(s, 1e-4), zeta_loc - 1e-4);
end

W        = best_W;
r        = best_r;
Rnew     = best_Rnew;
gaps_out = best_gaps;

if isempty(best_gaps)
    gaps_out = ones(8, 1);
end

if ~converged && ~quiet
    fprintf('    NOT converged after %d evals (best mkt_obj=%.3f)\n', n_evals, best_obj);
end
end

function g = eval_gaps(W, r, Rnew, struct_exp)
global gaps_calibration_last calib_mode hybrid_inner_clearing_skip

x_full = [W; r; Rnew; struct_exp];

hybrid_inner_clearing_skip = true;
calib_mode_old             = calib_mode;
calib_mode                 = true;
try
    evalc('evaluate_economy(x_full);');
catch ME_ic
    calib_mode                 = calib_mode_old;
    hybrid_inner_clearing_skip = false;
    rethrow(ME_ic);
end
calib_mode                 = calib_mode_old;
hybrid_inner_clearing_skip = false;

g = gaps_calibration_last;
end