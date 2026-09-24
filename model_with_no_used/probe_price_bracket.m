function [W_lo, W_hi, r_lo, r_hi, probe] = probe_price_bracket(theta_in, Wgrid, rgrid, verbose)
global eta_span eta_p omega kappa sigma_y exit_prob lambda_u beta ...
    gaps_calibration_last LAST_MOMENTS

if nargin < 4 || isempty(verbose), verbose = true; end

nW  = numel(Wgrid);
nR  = numel(rgrid);
gL  = nan(nW, nR);
gKn = nan(nW, nR);
gKK = nan(nW, nR);
ent = nan(nW, nR);

r_ceiling = 1/beta - 1;
if max(rgrid) >= r_ceiling
    warning('probe_price_bracket: rgrid reaches %.4f >= 1/beta-1 = %.4f (no stationary eq above it).', ...
        max(rgrid), r_ceiling);
end

if verbose
    fprintf('\n  [probe] coarse (W,r) grid at theta = %.4g   (%d x %d = %d evals)\n', ...
        theta_in, nW, nR, nW*nR);
    fprintf('          Rnew pinned in closed form -- this is a 2-price problem.\n\n');
end

jmid0 = max(1, round(nR/2));
Rn0   = rnew_closed_form(rgrid(jmid0));
x0    = [Wgrid(1); rgrid(jmid0); Rn0; 0; theta_in; eta_span; eta_p; omega; kappa; ...
    sigma_y;
    exit_prob;
    lambda_u];
evaluate_economy(x0);
g0 = gaps_calibration_last;
if verbose
    fprintf('    [pre-check] W=%.3f r=%.4f | L=%+8.1f%% Kn=%+8.1f%% ent=%.3f\n\n', ...
        Wgrid(1), rgrid(jmid0), 100*g0(2), 100*g0(1), LAST_MOMENTS.ent);
end
if g0(2) < 0
    error(['probe_price_bracket: the labour gap is ALREADY NEGATIVE (%+.1f%%) at the LOWEST ' ...
        'wage in the grid, W=%.3f. The root is BELOW the grid, so every point in it ' ...
        'would be on the wrong side of the equilibrium. Lower Wgrid -- start it just ' ...
        'above the CALIBRATED wage, not far above it. Labour demand is extremely ' ...
        'elastic in W here (occupational choice is a knife-edge under free credit), ' ...
        'so a huge labour gap does NOT imply a large wage adjustment.'], ...
        100*g0(2), Wgrid(1));
end

for i = 1:nW
    for j = 1:nR
        WW = Wgrid(i);
        rr = rgrid(j);
        Rn = rnew_closed_form(rr);
        x_full = [WW; rr; Rn; 0; theta_in; eta_span; eta_p; omega; kappa; ...
            sigma_y;
            exit_prob;
            lambda_u];
        evaluate_economy(x_full);
        g        = gaps_calibration_last;
        gKn(i,j) = g(1);
        gL(i,j)  = g(2);
        gKK(i,j) = g(3);
        ent(i,j) = LAST_MOMENTS.ent;
        if verbose
            fprintf('    W=%.3f r=%.4f Rnew=%.4f | L=%+8.1f%% Kn=%+8.1f%% KnKu=%+7.1f%% ent=%.3f\n', ...
                WW, rr, Rn, 100*gL(i,j), 100*gKn(i,j), 100*gKK(i,j), ent(i,j));
        end
    end
end

if verbose
    fprintf('\n  [probe] SIGN of the LABOUR gap  (+ = excess demand, wage too low)\n');
    fprintf('          rows = W, cols = r\n          %8s', '');
    fprintf(' %7.4f', rgrid);  fprintf('\n');
    for i = 1:nW
        fprintf('          W=%.3f ', Wgrid(i));
        for j = 1:nR
            if gL(i,j) > 0, fprintf('  %6s', '+'); else, fprintf('  %6s', '-'); end
        end
        fprintf('\n');
    end
    fprintf('\n  [probe] SIGN of the NEW-CAPITAL gap  (+ = excess demand, r too low)\n');
    fprintf('          %8s', '');
    fprintf(' %7.4f', rgrid);  fprintf('\n');
    for i = 1:nW
        fprintf('          W=%.3f ', Wgrid(i));
        for j = 1:nR
            if gKn(i,j) > 0, fprintf('  %6s', '+'); else, fprintf('  %6s', '-'); end
        end
        fprintf('\n');
    end
end

jmid  = max(1, round(nR/2));
W_lo  = Wgrid(1);
W_hi  = Wgrid(end);
iflip = find(gL(:,jmid) > 0, 1, 'last');
if isempty(iflip)
    warning(['probe_price_bracket: labour gap is NEGATIVE even at the LOWEST wage %.3f. ' ...
        'W_lo is too high -- extend Wgrid downward.'], Wgrid(1));
elseif iflip == nW
    warning(['probe_price_bracket: labour gap is STILL POSITIVE at the HIGHEST wage %.3f ' ...
        '(gap = %+.1f%%). The frictionless economy wants more labour than exists ' ...
        'at any wage in the grid -- extend Wgrid upward.'], Wgrid(end), 100*gL(nW,jmid));
else
    W_lo = Wgrid(iflip);
    W_hi = Wgrid(iflip+1);
end

irow = max(1, min(nW, iflip));
if isempty(iflip), irow = max(1, round(nW/2)); end
r_lo  = rgrid(1);
r_hi  = rgrid(end);
jflip = find(gKn(irow,:) > 0, 1, 'last');
if isempty(jflip)
    warning(['probe_price_bracket: capital gap is NEGATIVE even at the LOWEST r %.4f. ' ...
        'Extend rgrid downward.'], rgrid(1));
elseif jflip == nR
    warning(['probe_price_bracket: capital gap is STILL POSITIVE at the HIGHEST r %.4f ' ...
        '(gap = %+.1f%%). Capital demand exceeds savings at every admissible rate; ' ...
        'r cannot go above 1/beta-1 = %.4f, so there may be NO stationary ' ...
        'equilibrium.'], rgrid(end), 100*gKn(irow,nR), r_ceiling);
else
    r_lo = rgrid(jflip);
    r_hi = rgrid(jflip+1);
end

if verbose
    fprintf('\n  [probe] brackets: W in [%.3f %.3f],  r in [%.4f %.4f]\n', W_lo, W_hi, r_lo, r_hi);
    fprintf('          (pass these straight to clear_frictionless_bisect)\n\n');
end

probe = struct('Wgrid',Wgrid,'rgrid',rgrid,'gL',gL,'gKn',gKn,'gKK',gKK,'ent',ent, ...
    'theta',theta_in);
end
