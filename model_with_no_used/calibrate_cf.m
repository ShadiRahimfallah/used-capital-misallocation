function calibrate_cf(beta_free, n_refine, trust_frac, tol)
if nargin < 1 || isempty(beta_free),   beta_free  = false; end
if nargin < 2 || isempty(n_refine),    n_refine   = 4;     end
if nargin < 3 || isempty(trust_frac),  trust_frac = 0.25;  end
if nargin < 4 || isempty(tol),         tol        = 0.004; end

tag = 'bfix'; if beta_free, tag = 'bfree'; end
here = fileparts(mfilename('fullpath'));
addpath(here, '-end');

diary(sprintf('calibrate_cf_%s.log', tag)); diary on;
fprintf('\n==== calibrate_cf (%s)  %s ====\n\n', tag, datestr(now));

BETA_FIXED = 0.857687323825;

if beta_free
    PAR = {'theta','eta_p','omega','sigma_y','beta'};
    LB  = [1.40, 5.00, 0.72, 0.10, 0.8220];
    UB  = [3.40, 7.40, 0.95, 0.50, 0.8580];
    MOM = {'DtoY','wGini','exit','top10','r'};
    TGT = [0.6000, 0.3720, 0.09782, 0.6321, 0.0310];
else
    PAR = {'theta','eta_p','omega','sigma_y'};
    LB  = [1.40, 5.00, 0.72, 0.10];
    UB  = [3.40, 7.40, 0.95, 0.50];
    MOM = {'DtoY','wGini','exit','top10'};
    TGT = [0.6000, 0.3720, 0.09782, 0.6321];
end
np  = numel(PAR);
nm  = numel(MOM);
WID = UB - LB;

p0 = [2.14252333898, 5.94875930477, 0.86046679022, 0.244540551434];
if beta_free
    p0(5) = BETA_FIXED;

    if exist('calibrate_cf_bfix_result.mat','file')
        Q       = load('calibrate_cf_bfix_result.mat');
        p0(1:4) = Q.p_best(1:4);
        fprintf('seeded from the beta-fixed answer: theta=%.4f eta_p=%.4f omega=%.4f sigma_y=%.4f\n', p0(1:4));
    end
end

score    = @(m) sqrt(mean(((m - TGT)./abs(TGT)).^2));
seedfile = sprintf('surface_seed_cf_%s.mat', tag);

    function m = evalp(p, lbl)
        o = struct('theta',p(1),'eta_p',p(2),'omega',p(3),'sigma_y',p(4), ...
            'F',0.0,'rho_y',0.95,'do_clear',true,'fast',false, ...
            'W0',seedW,'r0',seedR,'Rnew0',seedRN,'tag',lbl);
        if beta_free, o.beta = p(5); else, o.beta = BETA_FIXED; end
        r_ = solve_equilibrium(o);
        if (isfield(r_,'failed') && r_.failed) || ~r_.converged
            fprintf('   !! %s did not clear -- dropped\n', lbl);  m = [];  return
        end
        seedW = r_.W; seedR = r_.r; seedRN = r_.Rnew;
        m     = [r_.DtoY, r_.wGini, r_.exit, r_.top10];
        if beta_free, m(5) = r_.r; end
        fprintf('   %s: ', lbl);
        for q = 1:nm, fprintf('%s=%.4f ', MOM{q}, m(q)); end
        fprintf('| RMS=%.4f\n', score(m));
    end

seedW = 0.550541310384; seedR = 0.0313782468366; seedRN = 0.0913782468366;
for sf = {'nightcf_cf_base.mat','run_cf_pure_result.mat'}
    if exist(sf{1},'file')
        z  = load(sf{1});
        fn = fieldnames(z);
        for q = 1:numel(fn)
            v = z.(fn{q});
            if isstruct(v) && isfield(v,'W') && isfield(v,'converged') && v.converged
                seedW  = v.W;
                seedR  = v.r;
                seedRN = v.Rnew;
                fprintf('seeded prices from %s: W=%.6f r=%.6f Rnew=%.6f\n\n', sf{1}, seedW, seedR, seedRN);
                break
            end
        end
        break
    end
end

if exist(seedfile, 'file')
    Sd = load(seedfile);
    X  = Sd.X;
    Y  = Sd.Y;
    fprintf('loaded %d cached design points from %s\n', size(X,1), seedfile);
else
    fprintf('building one-at-a-time design (%d evaluations)...\n', np+1);
    X  = [];
    Y  = [];
    m0 = evalp(p0, 'cf_base');
    if isempty(m0), error('calibrate_cf: base point did not clear -- fix the seeds first.'); end
    X(end+1,:) = p0;
    Y(end+1,:) = m0;
    for j = 1:np
        pj    = p0;
        pj(j) = min(max(p0(j) + 0.05*WID(j), LB(j)), UB(j));
        mj    = evalp(pj, sprintf('cf_oat_%s', PAR{j}));
        if ~isempty(mj), X(end+1,:) = pj; Y(end+1,:) = mj; end
    end
    save(seedfile,'X','Y','PAR','MOM');
end

s_all = arrayfun(@(i) score(Y(i,:)), 1:size(Y,1))';
[best_s, ib] = min(s_all);
p_best = X(ib,:);
m_best = Y(ib,:);
fprintf('\nbest design point: RMS=%.4f\n', best_s);

tf = trust_frac;
for it = 1:n_refine
    A = [ones(size(X,1),1), X];
    B = A \ Y;
    J = B(2:end,:)';
    if rank(J) < min(nm,np)
        fprintf('[iter %d] Jacobian rank-deficient -- stopping.\n', it); break
    end
    dp    = (J \ (-(m_best - TGT))')';
    sc    = min(1, min((tf*WID) ./ max(abs(dp),1e-12)));
    p_new = min(max(p_best + dp*sc, LB), UB);

    fprintf('\n---- refine %d (trust=%.3f) ----\n', it, tf);
    for j = 1:np, fprintf('   %-8s %8.4f -> %8.4f\n', PAR{j}, p_best(j), p_new(j)); end

    m_new = evalp(p_new, sprintf('cf_refine%d', it));
    if isempty(m_new), tf = tf/2; continue; end
    X(end+1,:) = p_new; Y(end+1,:) = m_new;
    s_new      = score(m_new);
    if s_new < best_s
        best_s = s_new;
        p_best = p_new;
        m_best = m_new;
        fprintf('   ACCEPTED (RMS -> %.4f)\n', best_s);
        tf = min(trust_frac, tf*1.5);
    else
        fprintf('   rejected (RMS %.4f) -- halving trust region\n', s_new);
        tf = tf/2;
    end
    save(seedfile,'X','Y','PAR','MOM');
    if best_s < tol
        fprintf('\n*** converged: RMS %.4f is below the tolerance %.4f ***\n', best_s, tol);
        break
    end
end

fprintf('\n==== calibrate_cf FINAL (%s) ====\n', tag);
fprintf('RMS relative gap = %.4f\n', best_s);
for j = 1:np, fprintf('   %-8s = %.4f\n', PAR{j}, p_best(j)); end
if ~beta_free, fprintf('   %-8s = %.4f  (FIXED at the benchmark value)\n','beta',BETA_FIXED); end
fprintf('\n   %-8s %10s %10s %9s\n','moment','model','target','gap%');
for i = 1:nm
    fprintf('   %-8s %10.4f %10.4f %+8.1f%%\n', MOM{i}, m_best(i), TGT(i), 100*(m_best(i)-TGT(i))/abs(TGT(i)));
end
save(sprintf('calibrate_cf_%s_result.mat',tag),'p_best','m_best','best_s','X','Y','PAR','MOM','TGT','beta_free');
fprintf('\n==== calibrate_cf DONE %s ====\n', datestr(now));
diary off;
end
