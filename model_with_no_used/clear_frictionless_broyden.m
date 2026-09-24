function [Wc, rc, Rnewc, converged, gaps_out, n_evals] = clear_frictionless_broyden( ...
    theta_in, W0, r0, J0, tol, max_evals, box, verbose)
global beta eta_span eta_p omega kappa sigma_y exit_prob lambda_u ...
    gaps_calibration_last LAST_MOMENTS zeta phi ce cw cxe cxw cse csw

HAS_USED = ~(isempty(zeta) || zeta <= 0 || phi >= 1 - 1e-12);

if nargin < 8 || isempty(verbose), verbose = true;  end
if nargin < 7 || isempty(box),     box = [0.40 0.80 0.020 0.154]; end

USE_FD = false;
if nargin >= 4 && ischar(J0) && strcmpi(J0, 'fd')
    USE_FD = true;
    J0     = [];
end
if isempty(J0)
    J0 = [ -9.0,  -6.1 ;
        -3.5, -22.6 ];
end

FD_DW = 0.004;
FD_DR = 0.002;

tolv = tol(:);
if isscalar(tolv), tolv = [tolv; tolv]; end
cvg = @(FF) (abs(FF(1)) < tolv(1)) && (abs(FF(2)) < tolv(2));

r_ceiling = 1/beta - 1;
box(4)    = min(box(4), 0.99*r_ceiling);

MAX_DW_FRAC = 0.02;
MAX_DR      = 0.005;
MAX_HALVE   = 3;

n_evals = 0;

    function [Fv, Rn, entv] = evalF(WW, rr)
        Rn = rnew_closed_form(rr);
        x_full = [WW; rr; Rn; 0; theta_in; eta_span; eta_p; omega; kappa; ...
            sigma_y;
            exit_prob;
            lambda_u];
        evaluate_economy(x_full);
        g       = gaps_calibration_last;
        Fv      = [g(2); g(1)];
        entv    = LAST_MOMENTS.ent;
        n_evals = n_evals + 1;
        if verbose
            if HAS_USED
                fprintf('    [broyden %2d] W=%.4f r=%.4f Rnew=%.4f | L=%+8.2f%% Kn=%+8.2f%% KnKu=%+6.2f%% ent=%.3f  ||F||=%.4f\n', ...
                    n_evals, WW, rr, Rn, 100*g(2), 100*g(1), 100*g(3), entv, norm(Fv));
            else
                fprintf('    [broyden %2d] W=%.4f r=%.4f Rnew=%.4f | L=%+8.2f%% Kn=%+8.2f%% ent=%.3f  ||F||=%.4f\n', ...
                    n_evals, WW, rr, Rn, 100*g(2), 100*g(1), entv, norm(Fv));
            end
        end
    end

    function xc = clamp_box(xx)
        xc    = xx;
        xc(1) = min(max(xc(1), box(1)), box(2));
        xc(2) = min(max(xc(2), box(3)), box(4));
    end

if verbose
    fprintf('\n  [clear_broyden] 2-D Broyden at theta = %.4g -- W and r move TOGETHER\n', theta_in);
    fprintf('    seed: W=%.4f r=%.4f   tol: Labor=%.2f%% Knew=%.2f%%   budget=%d evals\n', ...
        W0, r0, 100*tolv(1), 100*tolv(2), max_evals);
    fprintf('    J0 = [%6.2f %6.2f ; %6.2f %6.2f]   (measured; Broyden refines it)\n', ...
        J0(1,1), J0(1,2), J0(2,1), J0(2,2));
    fprintf('    Rnew pinned in closed form. Steps clipped to %.0f%% in W, %.3f in r.\n\n', ...
        100*MAX_DW_FRAC, MAX_DR);
end

x = clamp_box([W0; r0]);
[F, ~, ent] = evalF(x(1), x(2));
J = J0;

best_x    = x;
best_F    = F;
best_norm = norm(F, inf);

VF_SAVE = struct('ce',ce, 'cw',cw, 'cxe',cxe, 'cxw',cxw, 'cse',cse, 'csw',csw);

    function restore_vf()
        ce  = VF_SAVE.ce;
        cw  = VF_SAVE.cw;
        cxe = VF_SAVE.cxe;
        cxw = VF_SAVE.cxw;
        cse = VF_SAVE.cse;
        csw = VF_SAVE.csw;
    end

    function snapshot_vf()
        VF_SAVE = struct('ce',ce, 'cw',cw, 'cxe',cxe, 'cxw',cxw, 'cse',cse, 'csw',csw);
    end

converged      = cvg(F);
n_collapse_run = 0;

if USE_FD && ~converged
    if verbose
        fprintf('\n    [fd] measuring the Jacobian locally (2 evals). Forward steps only --\n');
        fprintf('         stepping r DOWN is what falls off the cliff.\n');
    end
    xw = clamp_box([x(1) + FD_DW; x(2)]);
    [Fw, ~, ent_w] = evalF(xw(1), xw(2));
    xr = clamp_box([x(1); x(2) + FD_DR]);
    [Fr, ~, ent_r] = evalF(xr(1), xr(2));

    if ent_w < 0.02 || ent_w > 0.95 || ent_r < 0.02 || ent_r > 0.95
        warning(['clear_frictionless_broyden: a finite-difference probe COLLAPSED the ' ...
            'entrepreneur share (ent_W=%.3f, ent_r=%.3f). The seed sits on the edge ' ...
            'of the basin. Keeping the default Jacobian; expect a fallback.'], ent_w, ent_r);
    else
        J = [ (Fw - F)/(xw(1) - x(1)), (Fr - F)/(xr(2) - x(2)) ];
        if abs(det(J)) < 1e-6
            warning('clear_frictionless_broyden: measured Jacobian is near-singular. Using the default.');
            J = J0;
        elseif verbose
            fprintf('    [fd] measured J = [%7.2f %7.2f ; %7.2f %7.2f]   det = %.1f\n', ...
                J(1,1), J(1,2), J(2,1), J(2,2), det(J));
            if sign(J(2,2)) ~= sign(J0(2,2))
                fprintf('    [fd] NOTE: dgKn/dr has the OPPOSITE sign to the frictionless case.\n');
                fprintf('         That is the collateral channel -- a lower r means less saving,\n');
                fprintf('         less wealth, less collateral, so SMALLER firms and LESS capital\n');
                fprintf('         demand. Assuming the frictionless sign here is what fell off the cliff.\n');
            end
        end
    end
end

while ~converged && n_evals < max_evals
    dx = -J \ F;
    if any(~isfinite(dx))
        warning('clear_frictionless_broyden: Jacobian went singular. Resetting to J0.');
        J  = J0;
        dx = -J \ F;
    end

    dx(1) = sign(dx(1)) * min(abs(dx(1)), MAX_DW_FRAC * x(1));
    dx(2) = sign(dx(2)) * min(abs(dx(2)), MAX_DR);

    accepted = false;
    for h = 0:MAX_HALVE
        step  = dx / 2^h;
        x_new = clamp_box(x + step);
        if norm(x_new - x) < 1e-9, break; end

        [F_new, ~, ent_new] = evalF(x_new(1), x_new(2));

        if ent_new < 0.02 || ent_new > 0.95
            if verbose
                fprintf('      -> rejected: entrepreneur share collapsed to %.3f (basin jump). Halving.\n', ent_new);
            end
            restore_vf();
            continue
        end
        if norm(F_new) < norm(F) || cvg(F_new)
            accepted = true;
            break
        end
        if verbose
            fprintf('      -> rejected: ||F|| rose %.4f -> %.4f. Halving.\n', norm(F), norm(F_new));
        end
        restore_vf();
        if n_evals >= max_evals, break; end
    end

    if ~accepted
        n_collapse_run = n_collapse_run + 1;
        restore_vf();

        if n_collapse_run >= 2 || n_evals >= max_evals
            warning(['clear_frictionless_broyden: boxed in after %d failed line search(es). ' ...
                'Every admissible step either collapses the entrepreneur share or ' ...
                'raises ||F||. Returning the best point: Labor %+.2f%%, Knew %+.2f%%. ' ...
                'If those gaps are large this parameter point may have NO equilibrium ' ...
                'on the healthy branch.'], ...
                n_collapse_run, 100*best_F(1), 100*best_F(2));
            break
        end

        J = J0;
        if verbose
            fprintf('      -> line search failed. Resetting J and retrying with a smaller step.\n');
        end
        MAX_DW_FRAC = MAX_DW_FRAC / 3;
        MAX_DR      = MAX_DR / 3;
        continue
    end

    s  = x_new - x;
    yv = F_new - F;
    if s' * s > 1e-14
        J = J + ((yv - J*s) * s') / (s' * s);
    end

    x = x_new;  F = F_new;  ent = ent_new;
    snapshot_vf();
    n_collapse_run = 0;
    if norm(F, inf) < best_norm
        best_x    = x;
        best_F    = F;
        best_norm = norm(F, inf);
    end
    converged = cvg(F);
end

if ~converged
    x = best_x;
    F = best_F;

    evalF(x(1), x(2));
end

Wc       = x(1);
rc       = x(2);
Rnewc    = rnew_closed_form(rc);
gaps_out = [F(2); F(1); NaN];

g_last   = gaps_calibration_last;
gaps_out = [g_last(1); g_last(2); g_last(3)];

if verbose
    if converged
        status = 'CONVERGED';
    else
        status = 'NOT CONVERGED';
    end
    fprintf('\n  [clear_broyden] %s after %d evals: W=%.4f r=%.4f Rnew=%.4f\n', ...
        status, n_evals, Wc, rc, Rnewc);
    if HAS_USED
        fprintf('                  gaps: Knew %+.2f%%  Labor %+.2f%%  KnKu %+.2f%%\n', ...
            100*gaps_out(1), 100*gaps_out(2), 100*gaps_out(3));
    else
        fprintf('                  gaps: Knew %+.2f%%  Labor %+.2f%%   (no used-capital market)\n', ...
            100*gaps_out(1), 100*gaps_out(2));
    end
    fprintf('                  final J = [%6.2f %6.2f ; %6.2f %6.2f]\n', ...
        J(1,1), J(1,2), J(2,1), J(2,2));
    if abs(rc - box(4)) < 1e-4
        fprintf('  [WARN] r pinned at the 1/beta-1 wall -- the capital market did NOT clear.\n');
    end
end
end
