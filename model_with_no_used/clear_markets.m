function [W, r, Rnew, converged, gaps_out] = clear_markets(struct_params, W0, r0, Rnew0, lb_p, ub_p)
global delta nested_clearing_quiet nested_clearing_max_iter ...
    nested_clearing_tol nested_clearing_tol_Knew

if isempty(delta), delta_local = 0.06; else, delta_local = delta; end

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

tol_kn = tol;
if ~isempty(nested_clearing_tol_Knew) && nested_clearing_tol_Knew > 0
    tol_kn = nested_clearing_tol_Knew;
end

lb         = lb_p(:);
ub         = ub_p(:);
step_limit = 0.15;
fd_rel     = 0.02;
max_bt     = 5;

n_evals = 0;

x = proj([W0; r0], lb, ub);
[F, g, ok] = resid_c(x, struct_exp);
if ~ok
    x = proj([mean([lb(1) ub(1)]); mean([lb(2) ub(2)])], lb, ub);
    [F, g, ok] = resid_c(x, struct_exp);
    if ~ok
        W         = x(1);
        r         = x(2);
        Rnew      = r + delta_local;
        converged = false;
        gaps_out  = ones(6,1);
        if ~quiet, fprintf('  [clear] no admissible starting point -- returning penalty gaps\n'); end
        return;
    end
end

best_norm = max(abs(F));
best_x    = x;
best_g    = g;

cvg       = @(FF) (abs(FF(1)) < tol) && (abs(FF(2)) < tol_kn);
converged = cvg(F);

if ~quiet
    fprintf('\n  [clear] TWO-MARKET Broyden (W: labor, r: new capital; Rnew=r+delta pinned)\n');
    fprintf('  [clear] budget %d evals, tol: Labor=%.1f%%  Knew=%.1f%% (cliff band)\n', ...
        max_evals, 100*tol, 100*tol_kn);
    fprintf('  [clear] start  W=%.4f r=%.4f Rnew=%.4f | L=%+.1f%% Kn=%+.1f%%\n', ...
        x(1), x(2), x(2)+delta_local, 100*F(1), 100*F(2));
end

J = -eye(2);
if ~converged
    J = fd_jacobian_c(x, F, struct_exp, lb, ub, fd_rel);
end

while ~converged && n_evals < max_evals
    if abs(det(J)) < 1e-9, J = J + 1e-2*eye(2); end

    dx = -(J \ F);
    if any(~isfinite(dx)), dx = -0.1 * F; end
    for j = 1:2
        cap = step_limit * max(abs(x(j)), 1e-3);
        if abs(dx(j)) > cap, dx(j) = sign(dx(j)) * cap; end
    end

    lam      = 1.0;
    accepted = false;
    x_new    = x;
    F_new    = F;
    g_new    = g;
    for bt = 0:max_bt
        if n_evals >= max_evals, break; end
        xt = proj(x + lam*dx, lb, ub);
        [Ft, gt, okt] = resid_c(xt, struct_exp);
        if okt && all(isfinite(Ft)) && max(abs(Ft)) < max(abs(F)) * (1 - 1e-4)
            x_new    = xt;
            F_new    = Ft;
            g_new    = gt;
            accepted = true;
            break;
        end
        lam = 0.5 * lam;
    end

    if ~accepted
        if n_evals >= max_evals, break; end
        J = fd_jacobian_c(x, F, struct_exp, lb, ub, fd_rel);
        if abs(det(J)) < 1e-9, J = J + 1e-2*eye(2); end
        dx = -(J \ F);
        for j = 1:2
            cap = step_limit * max(abs(x(j)), 1e-3);
            if abs(dx(j)) > cap, dx(j) = sign(dx(j)) * cap; end
        end
        xt = proj(x + 0.5*dx, lb, ub);
        [Ft, gt, okt] = resid_c(xt, struct_exp);
        if okt && all(isfinite(Ft))
            x_new = xt;
            F_new = Ft;
            g_new = gt;
        else
            if ~quiet, fprintf('  [clear] no admissible step -- stop (best max|F|=%.4f)\n', best_norm); end
            break;
        end
    end

    s = x_new - x;  yv = F_new - F;  sts = s.'*s;
    if sts > 1e-14
        J = J + ((yv - J*s) * s.') / sts;
    end

    x = x_new;
    F = F_new;
    g = g_new;
    if max(abs(F)) < best_norm
        best_norm = max(abs(F));
        best_x    = x;
        best_g    = g;
    end
    converged = cvg(F);

    if ~quiet
        msg = '';  if converged, msg = '   CONVERGED'; end
        fprintf('  [clear] W=%.4f r=%.4f Rnew=%.4f | L=%+.1f%% Kn=%+.1f%% | max|F|=%.4f (evals %d/%d)%s\n', ...
            x(1), x(2), x(2)+delta_local, 100*F(1), 100*F(2), max(abs(F)), n_evals, max_evals, msg);
    end
end

W    = best_x(1);
r    = best_x(2);
Rnew = r + delta_local;
if isempty(best_g) || numel(best_g) < 6
    gaps_out  = ones(6,1);
    converged = false;
else
    gaps_out  = best_g;

    converged = (abs(best_g(1)) < tol_kn) && (abs(best_g(2)) < tol);
end

if ~quiet && ~converged
    fprintf('  [clear] NOT converged after %d evals (best max|F|=%.4f) W=%.4f r=%.4f\n', ...
        n_evals, best_norm, W, r);
end

    function xp = proj(xin, lbv, ubv)
        xp    = xin;
        xp(1) = min(max(xp(1), lbv(1)), ubv(1));
        xp(2) = min(max(xp(2), lbv(2)), ubv(2));
    end

    function [Fv, gv, okv] = resid_c(xin, sexp)
        global gaps_calibration_last calib_mode hybrid_inner_clearing_skip ...
            LAST_MOMENTS

        x_full = [xin(1); xin(2); xin(2) + delta_local; sexp];

        hisk                       = hybrid_inner_clearing_skip;
        cmo                        = calib_mode;
        hybrid_inner_clearing_skip = true;
        calib_mode                 = true;
        try
            evaluate_economy(x_full);
        catch
            hybrid_inner_clearing_skip = hisk;
            calib_mode                 = cmo;
            n_evals                    = n_evals + 1;
            Fv                         = [NaN; NaN];
            gv                         = [];
            okv                        = false;
            return;
        end
        hybrid_inner_clearing_skip = hisk;
        calib_mode                 = cmo;
        n_evals                    = n_evals + 1;

        gv = gaps_calibration_last;

        if isempty(gv) || numel(gv) < 6 || any(abs(gv(1:2)) > 3)
            Fv  = [NaN; NaN];
            gv  = [];
            okv = false;
            return;
        end

        if ~isempty(LAST_MOMENTS) && isstruct(LAST_MOMENTS) && isfield(LAST_MOMENTS,'ent')
            ent_v = LAST_MOMENTS.ent;
            if ~isfinite(ent_v) || ent_v < 0.03 || ent_v > 0.90
                Fv  = [NaN; NaN];
                gv  = [];
                okv = false;
                return;
            end
        end
        Fv  = [gv(2); gv(1)];
        okv = true;
    end

    function Jout = fd_jacobian_c(xin, Fin, sexp, lbv, ubv, fdrel)
        Jout = -eye(2);
        for j = 1:2
            if n_evals >= max_evals, return; end
            dxj   = fdrel * max(abs(xin(j)), 1e-3);
            xj    = xin;
            xj(j) = xj(j) + dxj;
            xj    = proj(xj, lbv, ubv);
            h     = xj(j) - xin(j);
            if abs(h) < 1e-10
                xj    = xin;
                xj(j) = xj(j) - dxj;
                xj    = proj(xj, lbv, ubv);
                h     = xj(j) - xin(j);
            end
            if abs(h) < 1e-10, continue; end
            [Fj, ~, okj] = resid_c(xj, sexp);
            if okj && all(isfinite(Fj))
                Jout(:, j) = (Fj - Fin) / h;
            end
        end
    end
end
