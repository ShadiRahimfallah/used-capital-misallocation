function out = solve_equilibrium_cached(opts)
if nargin < 1, opts = struct(); end
here  = fileparts(mfilename('fullpath'));
cfile = fullfile(here, 'equilibrium_cache.mat');

key = cache_key(opts);

keys = {};
results = {};
if exist(cfile, 'file')
    try
        S = load(cfile);
        if isfield(S, 'keys') && isfield(S, 'results')
            keys    = S.keys;
            results = S.results;
        end
    catch ME
        warning('solve_equilibrium_cached: cache unreadable (%s) -- starting a new one.', ME.message);
    end
end

hit = find(strcmp(keys, key), 1);
if ~isempty(hit)
    out = results{hit};
    fprintf('  [cache HIT ] %s\n', key);
    fprintf('               r=%.4f exit=%.4f wGini=%.4f top10=%.4f DtoY=%.4f su=%+.4f\n', ...
        gf(out,'r'), gf(out,'exit'), gf(out,'wGini'), gf(out,'top10'), gf(out,'DtoY'), gf(out,'beta_size'));
    return
end

fprintf('  [cache MISS] %s\n', key);
t0  = tic;
out = solve_equilibrium(opts);
fprintf('  [cache ADD ] evaluated in %.1f min\n', toc(t0)/60);

keys{end+1}    = key;
results{end+1} = out;
try
    save(cfile, 'keys', 'results');
catch ME
    warning('solve_equilibrium_cached: could not write the cache (%s).', ME.message);
end
end

function v = gf(s, f)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = NaN;
end
end
