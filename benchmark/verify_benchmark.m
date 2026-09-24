function verify_benchmark(fast_flag)
if nargin < 1 || isempty(fast_flag), fast_flag = false; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'benchmark'), '-end');
addpath(here, '-end');

diary(fullfile(here, 'verify_benchmark.log')); diary on;

PAR = {'theta','eta_p','omega','kappa','sigma_y','beta'};
MOM = {'r','exit','wGini','top10','DtoY','su'};
TGT = [0.0310, 0.09782, 0.3720, 0.6321, 0.6000, -0.0269];

p = read_published_parameters(fullfile(here, 'parameters_benchmark.m'), PAR);

fprintf('\n================================================================\n');
fprintf('  INDEPENDENT VERIFICATION OF THE PUBLISHED BENCHMARK\n');
fprintf('  %s\n', datestr(now));
fprintf('================================================================\n\n');
fprintf('  This script calls solve_equilibrium DIRECTLY. It does not read\n');
fprintf('  equilibrium_cache.mat and it does not read benchmark_result.mat\n');
fprintf('  for anything except the comparison printed at the end. Every moment\n');
fprintf('  below is recomputed from the parameters.\n\n');
fprintf('  grid = %s\n\n', ternary_str(fast_flag, 'k32 (fast, NOT the production grid)', 'k50 (production)'));
for j = 1:6
    fprintf('     %-8s = %.4f\n', PAR{j}, p(j));
end

opts = struct('theta',p(1),'eta_p',p(2),'omega',p(3),'kappa',p(4), ...
    'sigma_y',p(5),'beta',p(6),'F',0.0,'rho_y',0.95,'lambda_u',1.0, ...
    'do_clear',true,'fast',fast_flag,'tag','verify_benchmark');

fprintf('\n  solving (a cold solve on the production grid takes roughly 40 minutes)...\n\n');
t0 = tic;
o  = solve_equilibrium(opts);
el = toc(t0)/60;

m = [o.r, o.exit, o.wGini, o.top10, o.DtoY, o.beta_size];
rms = sqrt(mean(((m - TGT)./abs(TGT)).^2));

fprintf('\n================================================================\n');
fprintf('  RESULT  (solved in %.1f minutes, converged = %d)\n', el, o.converged);
fprintf('================================================================\n\n');
fprintf('  prices: W = %.6f   r = %.6f   Rnew = %.6f\n\n', o.W, o.r, o.Rnew);
fprintf('  %-8s %12s %12s %10s\n', 'moment', 'recomputed', 'target', 'gap%');
for i = 1:6
    fprintf('  %-8s %12.6f %12.6f %+9.1f%%\n', MOM{i}, m(i), TGT(i), ...
        100*(m(i)-TGT(i))/abs(TGT(i)));
end
fprintf('\n  RMS relative gap = %.4f\n', rms);

bfile = fullfile(here, 'benchmark_result.mat');
if exist(bfile, 'file')
    B = load(bfile);
    r0 = [B.out.r, B.out.exit, B.out.wGini, B.out.top10, B.out.DtoY, B.out.beta_size];
    fprintf('\n  COMPARISON AGAINST THE MOMENTS REPORTED IN THE PAPER\n\n');
    fprintf('  %-8s %12s %12s %14s\n', 'moment', 'recomputed', 'reported', 'difference');
    for i = 1:6
        fprintf('  %-8s %12.6f %12.6f %14.2e\n', MOM{i}, m(i), r0(i), m(i)-r0(i));
    end
    dmax = max(abs(m - r0));
    fprintf('\n  largest absolute difference = %.3e\n', dmax);
    if dmax < 1e-4
        fprintf('  VERDICT: the reported moments reproduce.\n');
    else
        fprintf('  VERDICT: DOES NOT REPRODUCE -- investigate before using the reported numbers.\n');
    end
    if fast_flag
        fprintf('\n  NOTE you ran the k32 grid. top10 is a tail statistic and is grid\n');
        fprintf('  sensitive, so a top10 difference here is expected and is not evidence\n');
        fprintf('  of a problem. Re-run with verify_benchmark(false) to compare properly.\n');
    end
end

save(fullfile(here, 'verify_benchmark_result.mat'), 'o', 'p', 'm', 'TGT', 'rms', 'fast_flag');
fprintf('\n  saved verify_benchmark_result.mat\n');

if o.converged
    cfile = fullfile(here, 'equilibrium_cache.mat');
    key   = cache_key(opts);
    keys = {}; results = {};
    if exist(cfile, 'file')
        try
            S = load(cfile);
            if isfield(S,'keys') && isfield(S,'results')
                keys = S.keys; results = S.results;
            end
        catch
        end
    end
    hit = find(strcmp(keys, key), 1);
    if isempty(hit)
        keys{end+1} = key; results{end+1} = o;
    else
        results{hit} = o;
    end
    try
        save(cfile, 'keys', 'results');
        fprintf('  wrote this solve into equilibrium_cache.mat, so recalibrate\n');
        fprintf('  will not spend another 40 minutes on the same point.\n');
    catch ME
        warning('verify_benchmark: could not update the cache (%s).', ME.message);
    end
end
fprintf('================================================================\n');
diary off;
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
