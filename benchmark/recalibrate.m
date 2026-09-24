function recalibrate(n_iter, fast_flag, trust_frac, tol)
if nargin < 1 || isempty(n_iter),     n_iter     = 6;     end
if nargin < 2 || isempty(fast_flag),  fast_flag  = false; end
if nargin < 3 || isempty(trust_frac), trust_frac = 0.10;  end
if nargin < 4 || isempty(tol),        tol        = 0.004; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'benchmark'), '-end');
addpath(here, '-end');

PAR = {'theta','eta_p','omega','kappa','sigma_y','beta'};
MOM = {'r','exit','wGini','top10','DtoY','su'};
TGT = [0.0310, 0.09782, 0.3720, 0.6321, 0.6000, -0.0269];

p_final = read_published_parameters(fullfile(here, 'parameters_benchmark.m'), PAR);

vfile = fullfile(here, 'verify_benchmark_result.mat');
bfile = fullfile(here, 'benchmark_result.mat');
anchor_src = '';
if exist(vfile, 'file')
    V = load(vfile);
    if isfield(V,'fast_flag') && ~V.fast_flag && isfield(V,'m')
        m_final    = V.m(:)';
        p_final    = V.p(:)';
        anchor_src = 'your own verify_benchmark run (k50)';
    end
end
if isempty(anchor_src)
    if ~exist(bfile, 'file')
        error('recalibrate: missing %s -- cannot anchor on the published point.', bfile);
    end
    B = load(bfile);
    o = B.out;
    m_final    = [o.r, o.exit, o.wGini, o.top10, o.DtoY, o.beta_size];
    p_final    = [o.theta, o.eta_p, o.omega, o.kappa, o.sigma_y, o.beta];
    anchor_src = 'benchmark_result.mat (run verify_benchmark first if you want to anchor on your own solve)';
end

score = sqrt(mean(((m_final - TGT)./abs(TGT)).^2));

fprintf('\n================================================================\n');
fprintf('  RECALIBRATION -- parameter search started at the published point\n');
fprintf('================================================================\n\n');
fprintf('  grid              : %s\n', ternary_str(fast_flag, 'k32 (fast)', 'k50 (production)'));
fprintf('  Newton iterations : %d\n', n_iter);
fprintf('  trust region      : %.3f of the box width per step\n', trust_frac);
fprintf('  stop when RMS <   : %.4f\n', tol);
fprintf('  evaluation cache  : %s\n', fullfile(here, 'equilibrium_cache.mat'));
fprintf('\n  starting parameters (published benchmark)\n');
for j = 1:6
    fprintf('     %-8s = %.4f\n', PAR{j}, p_final(j));
end
fprintf('\n  anchor source     : %s\n', anchor_src);
fprintf('\n  starting moments (read, not recomputed)\n');
fprintf('     %-8s %10s %10s %9s\n', 'moment', 'model', 'target', 'gap%');
for i = 1:6
    fprintf('     %-8s %10.4f %10.4f %+8.1f%%\n', MOM{i}, m_final(i), TGT(i), ...
        100*(m_final(i)-TGT(i))/abs(TGT(i)));
end
fprintf('\n  RMS relative gap at the starting point = %.4f   <== the number to beat\n', score);
fprintf('\n  The anchor costs no solve. Any Newton step that revisits a cached\n');
fprintf('  parameter vector returns instantly; only genuinely new points are\n');
fprintf('  solved, and each one is added to the cache as it completes.\n');
fprintf('\n  A step is ACCEPTED only if it lowers the RMS, so this search cannot\n');
fprintf('  return a worse point than the one it started from.\n\n');

calibrate_newton(n_iter, trust_frac, 0.60, fast_flag, p_final, m_final, tol);

fprintf('\n================================================================\n');
fprintf('  DONE. Result saved by calibrate_newton to\n');
fprintf('    calibrate_newton_result_%s.mat\n', ternary_str(fast_flag,'k32','k50'));
fprintf('  Cached evaluations persist in equilibrium_cache.mat, so re-running\n');
fprintf('  this script repeats no work.\n');
fprintf('================================================================\n');
end

function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end

function p = read_published_parameters(pfile, names)
if ~exist(pfile, 'file')
    error('read_published_parameters: %s not found.', pfile);
end
txt = fileread(pfile);
p   = nan(1, numel(names));
for j = 1:numel(names)
    tok = regexp(txt, ['^\s*' names{j} '\s*=\s*([^;%]+);'], ...
        'tokens', 'lineanchors', 'once');
    if isempty(tok)
        error('read_published_parameters: %s not found in %s.', names{j}, pfile);
    end
    p(j) = str2double(strtrim(tok{1}));
    if ~isfinite(p(j))
        error('read_published_parameters: %s is not a plain number in %s.', ...
            names{j}, pfile);
    end
end
end
