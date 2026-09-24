function calibrate_newton(n_iter, trust_frac, kernel_bw, fast_flag, p_start, m_start, tol)
if nargin < 1 || isempty(n_iter),     n_iter     = 5;    end
if nargin < 2 || isempty(trust_frac), trust_frac = 0.25; end
if nargin < 3 || isempty(kernel_bw),  kernel_bw  = 0.60; end
if nargin < 4 || isempty(fast_flag),  fast_flag  = true; end
if nargin < 5, p_start = []; end
if nargin < 6, m_start = []; end
if nargin < 7 || isempty(tol),        tol        = 0.02; end
gtag = 'k32'; if ~fast_flag, gtag = 'k50'; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'benchmark'), '-end');

diary(fullfile(here, 'calibrate_newton.log')); diary on;
fprintf('\n==== calibrate_newton  %s  (n_iter=%d trust=%.2f bw=%.2f grid=%s) ====\n\n', ...
    datestr(now), n_iter, trust_frac, kernel_bw, ternary_str(fast_flag,'k32','k50'));

PAR = {'theta','eta_p','omega','kappa','sigma_y','beta'};
MOM = {'r','exit','wGini','top10','DtoY','su'};
TGT = [0.0310, 0.09782, 0.3720, 0.6321, 0.6000, -0.0269];
LB  = [1.40, 5.00, 0.72, 0.01, 0.10, 0.8220];
UB  = [3.40, 7.40, 0.95, 0.20, 0.50, 0.8580];
WID = UB - LB;

BETA_SAFE_MAX = 0.8600;

S = load('calibrate_newton_seed.mat');
X = S.X;
Y = S.Y;
fprintf('loaded %d seed evaluations\n', size(X,1));

logfile = ['calibrate_newton_log_' gtag '.mat'];
if exist(logfile, 'file')
    try
        L = load(logfile);
        if isfield(L,'X') && isfield(L,'Y') && size(L.X,2) == size(X,2)
            n0 = size(X,1);
            for i = 1:size(L.X,1)
                if ~any(all(abs(X - L.X(i,:)) < 1e-12, 2))
                    X(end+1,:) = L.X(i,:);
                    Y(end+1,:) = L.Y(i,:);
                end
            end
            fprintf('merged %d further evaluations from %s (total %d)\n', ...
                size(X,1)-n0, logfile, size(X,1));
        end
    catch ME
        warning('calibrate_newton: could not merge %s (%s).', logfile, ME.message);
    end
end

score = @(m) sqrt(mean(((m - TGT)./abs(TGT)).^2));

s_all = arrayfun(@(i) score(Y(i,:)), 1:size(Y,1))';
[best_s, ib] = min(s_all);
p_best = X(ib,:);
m_best = Y(ib,:);
if ~isempty(p_start) && ~isempty(m_start)
    p_best = p_start(:)';  m_best = m_start(:)';  best_s = score(m_best);
    fprintf('ANCHORED on the supplied point (slopes from the seed set, levels from here).\n');
end
fprintf('best seed: RMS rel gap = %.4f at\n   ', best_s);
fprintf('%s=%.4f ', PAR{1}, p_best(1)); for j=2:6, fprintf('%s=%.4f ', PAR{j}, p_best(j)); end
fprintf('\n');

seedW  = 0.5125;
seedR  = m_best(1);
seedRN = 0.1109;

tf = trust_frac;
for it = 1:n_iter
    d  = sqrt(sum(((X - p_best)./WID).^2, 2));
    wt = exp(-0.5*(d/kernel_bw).^2);
    A  = [ones(size(X,1),1), X];
    Aw = A .* wt;
    Yw = Y .* wt;
    B  = Aw \ Yw;
    J  = B(2:end,:)';

    if rcond(J) < 1e-10
        fprintf('\n[iter %d] Jacobian near-singular (rcond=%.2e) -- stopping.\n', it, rcond(J));
        break
    end

    gap    = m_best - TGT;
    dp_raw = (J \ (-gap(:)))';
    cap    = tf * WID;
    ratios = cap ./ max(abs(dp_raw), 1e-12);
    [rmin, ibind] = min(ratios);
    scale  = min(1, rmin);
    dp     = dp_raw * scale;
    p_new  = min(max(p_best + dp, LB), UB);

    p_new(6) = min(p_new(6), BETA_SAFE_MAX);
    if p_best(6) + dp(6) > BETA_SAFE_MAX
        fprintf('   [beta guard] step wanted beta=%.4f, held at the safe ceiling %.4f\n', ...
            p_best(6) + dp(6), BETA_SAFE_MAX);
    end

    fprintf('\n---- iteration %d (trust=%.3f, step scaled by %.3f) ----\n', it, tf, scale);
    for j = 1:6
        fprintf('   %-8s %8.4f -> %8.4f   (%+.3f box widths)%s\n', ...
            PAR{j}, p_best(j), p_new(j), (p_new(j)-p_best(j))/WID(j), ...
            ternary_str(scale < 1 && j == ibind, '  <== binds the trust region', ''));
    end
    pred = m_best + (J*(p_new-p_best)')';
    fprintf('   predicted: ');
    for i = 1:6, fprintf('%s=%.4f ', MOM{i}, pred(i)); end
    fprintf('\n');

    opts = struct('theta',p_new(1),'eta_p',p_new(2),'omega',p_new(3), ...
        'kappa',p_new(4),'sigma_y',p_new(5),'beta',p_new(6), ...
        'F',0.0,'rho_y',0.95,'lambda_u',1.0,'do_clear',true,'fast',fast_flag, ...
        'W0',seedW,'r0',max(min(seedR,0.08),0.005),'Rnew0',seedRN, ...
        'tag',sprintf('search_%s_i%d',gtag,it));
    try
        o = solve_equilibrium_cached(opts);
    catch ME
        fprintf('   !! evaluation failed: %s -- halving trust region\n', ME.message);
        tf = tf/2;
        continue
    end
    if ~o.converged
        fprintf('   !! did not clear -- halving trust region, point discarded\n');
        tf = tf/2;
        continue
    end

    seedW      = o.W;  seedR = o.r;  seedRN = o.Rnew;
    m_new      = [o.r, o.exit, o.wGini, o.top10, o.DtoY, o.beta_size];
    X(end+1,:) = p_new;  Y(end+1,:) = m_new;
    s_new      = score(m_new);

    fprintf('   %-8s %10s %10s %10s %9s\n', 'moment', 'predicted', 'realised', 'target', 'gap%');
    for i = 1:6
        fprintf('   %-8s %10.4f %10.4f %10.4f %+8.1f%%\n', ...
            MOM{i}, pred(i), m_new(i), TGT(i), 100*(m_new(i)-TGT(i))/abs(TGT(i)));
    end
    if s_new < best_s, verdict = 'IMPROVED'; else, verdict = 'no improvement'; end
    fprintf('   RMS relative gap: %.4f -> %.4f  (%s)\n', best_s, s_new, verdict);

    if s_new < best_s
        best_s = s_new;
        p_best = p_new;
        m_best = m_new;
        fprintf('   ACCEPTED as new best.\n');
        tf = min(trust_frac, tf*1.5);
    else
        fprintf('   rejected (no improvement) -- halving trust region.\n');
        tf = tf/2;
    end

    save(['calibrate_newton_log_' gtag '.mat'], 'X', 'Y', 'PAR', 'MOM');
    if best_s < tol
        fprintf('\n*** converged: RMS %.4f is below the tolerance %.4f ***\n', best_s, tol);
        break
    end
end

fprintf('\n==== calibrate_newton FINAL ====\n');
fprintf('RMS relative gap = %.4f\n', best_s);
for j = 1:6, fprintf('   %-8s = %.4f\n', PAR{j}, p_best(j)); end
fprintf('\n   %-8s %10s %10s %9s\n', 'moment', 'model', 'target', 'gap%');
for i = 1:6
    fprintf('   %-8s %10.4f %10.4f %+8.1f%%\n', MOM{i}, m_best(i), TGT(i), ...
        100*(m_best(i)-TGT(i))/abs(TGT(i)));
end
save(['calibrate_newton_log_' gtag '.mat'], 'X', 'Y', 'PAR', 'MOM');
save(['calibrate_newton_result_' gtag '.mat'],'p_best','m_best','best_s','X','Y');
fprintf('\n==== calibrate_newton DONE %s ====\n', datestr(now));
diary off;
end

function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end
