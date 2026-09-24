function tables()
% tables -- every table the paper reports, printed in about a minute from
% results that are already stored in this package. Nothing is re-solved and
% no number in this file is typed in by hand.
%
%   >> tables
%
% Writes TABLES.txt next to this file.
%
%   Table 1  Used capital by firm age                       empirical
%   Table 2  Used capital by employment size                empirical
%   Table 3  Calibration: parameters, values and targets    model
%   Table 4  Used capital and capital misallocation         model
%   Table 5  Used capital and entry by productive poor      model
%   Table 6  The result under three designs without used K  model
%   Table 7  Sensitivity of the results to phi              model
%
% WHERE THE NUMBERS COME FROM
%   Tables 1 and 2 are echoed verbatim from empirical/motivation/
%   used_capital_table.txt, which the Python script in that folder wrote from
%   the Vietnam SME survey. That microdata is not included in this package,
%   so the committed results file is the audit trail.
%
%   Tables 3 to 7 are computed here, in this file, from the solved equilibria
%   stored in benchmark/, model_with_no_used/, cf_calib_B/ and cf_etap/. Every
%   cell is read from a .mat field or derived from fields by an expression you
%   can see in the source below. The appendix printed at the end lists each
%   source file with its size and modification date, and names the field
%   behind every row, so any cell can be traced back to the solve that
%   produced it.
%
%   To re-derive the stored equilibria from scratch instead of reading them,
%   run benchmark/verify_benchmark.m -- that takes about forty minutes.

HERE = fileparts(mfilename('fullpath'));
cd(HERE); addpath(HERE);

OUT = fullfile(HERE, 'TABLES.txt');
diary off;
try, if exist(OUT,'file'), delete(OUT); end, catch, end   %#ok<TRYNC>
diary(OUT); diary on;

L = @(c) fprintf('%s\n', repmat(c, 1, 78));

L('=');
fprintf('  FINANCIAL FRICTIONS, USED CAPITAL, AND MISALLOCATION\n');
fprintf('  Every table in the paper, printed from stored results\n');
fprintf('  generated %s\n', char(datetime('now')));
L('=');
fprintf('\n  This script re-solves nothing. It reads the saved equilibria and\n');
fprintf('  prints Tables 1 to 7. It runs in about a minute. The provenance\n');
fprintf('  appendix at the end names the file and field behind every number.\n');

TGT = paper_targets();

% ---- load every solved economy -------------------------------------------
SPEC = { ...
    'Benchmark',       'benchmark',          'benchmark_result.mat', 'analysis_used_capital_result.mat'
    'No used capital', 'model_with_no_used', 'cf_result.mat',        'analysis_cf_result.mat'
    'Tail pinned',     'cf_calib_B',         'cf_result.mat',        'analysis_cf_result.mat'
    'Tail diagnostic', 'cf_etap',            'cf_result.mat',        'analysis_cf_result.mat'
    };
[E, SRC] = load_all(SPEC);

ib = find(strcmp({E.name}, 'Benchmark'),       1);
ic = find(strcmp({E.name}, 'No used capital'), 1);
ip = find(strcmp({E.name}, 'Tail pinned'),     1);
id = find(strcmp({E.name}, 'Tail diagnostic'), 1);

if isempty(ib)
    diary off;
    error('tables:noBenchmark', ...
        ['The benchmark economy is missing. Expected benchmark/benchmark_result.mat ' ...
         'and benchmark/nofriction_result.mat. Nothing can be normalised without it.']);
end

missing = setdiff(SPEC(:,1)', {E.name});
if ~isempty(missing)
    fprintf('\n  ABSENT: %s -- the matching column or row is blank below.\n', ...
        strjoin(missing, ', '));
end


%% ------------------------------------------------------------------ T1, T2
echo_empirical(L, fullfile(HERE, 'empirical', 'motivation', 'used_capital_table.txt'));


%% ---------------------------------------------------------------- TABLE 3
head(L, 'TABLE 3.  Calibration: Parameters, Values, and Targets');

PB = read_params(fullfile(HERE,'benchmark','parameters_benchmark.m'));
PC = read_params(fullfile(HERE,'model_with_no_used','parameters_cf_calib.m'));

fprintf('\n  %-44s%11s%15s   %s\n', 'Parameter', 'Benchmark', 'Model without', 'Source or target');
fprintf('  %-44s%11s%15s\n', '', '', 'used capital');

fprintf('\n  Panel A: Parameters assigned from the literature\n');
arow('alpha    Labor share',              gp(PB,'alpha'),      gp(PC,'alpha'),      'Buera et al. (2011)');
arow('delta    Depreciation rate',        gp(PB,'delta'),      gp(PC,'delta'),      'Buera et al. (2011)');
arow('eta      Span of control',          gp(PB,'eta_span'),   gp(PC,'eta_span'),   'Buera et al. (2011)');
arow('sigma    Risk aversion',            gp(PB,'sigma_crra'), gp(PC,'sigma_crra'), 'Buera et al. (2011)');
arow('rho_y    Ability persistence',      gp(PB,'rho_y'),      gp(PC,'rho_y'),      'Cagetti and De Nardi (2006)');
arow('phi      CES weight on new capital',gp(PB,'phi'),        gp(PC,'phi'),        'Assigned; see text');
% gamma governs substitution between the two vintages. With phi = 1 there is
% only one vintage, so it is undefined there rather than equal to anything.
gC = gp(PC,'gamma'); if gp(PC,'phi') >= 1 - 1e-12, gC = NaN; end
arow('gamma    Substitution, new vs used',gp(PB,'gamma'),      gC,                  'Edgerton (2011)');

fprintf('\n  Panel B: Parameters jointly calibrated to data moments\n');
brow(E, ib, ic, 'theta',   'theta    Collateral multiplier',        sprintf('Debt-to-output ratio: %.3f',            TGT.DtoY));
brow(E, ib, ic, 'eta_p',   'eta_p    Productivity-tail parameter',  sprintf('Top-decile employment share: %.3f',     TGT.top10));
brow(E, ib, ic, 'omega',   'omega    Productivity persistence',     sprintf('Firm exit rate: %.3f',                  TGT.exit));
brow(E, ib, ic, 'sigma_y', 'sigma_y  Ability innovation dispersion',sprintf('Worker wage Gini: %.3f',                TGT.wGini));
brow(E, ib, ic, 'beta',    'beta     Discount factor',              sprintf('Real interest rate: %.3f',              TGT.r));

% zeta is implied by the stationary vintage-stock condition, not searched
% over: zeta = (delta + kappa) / (Kn/Ku). Where there is no used capital
% there is no transition, so it is zero rather than undefined.
zb = zeta_of(E, ib, TGT);   zc = zeta_of(E, ic, TGT);
prow3('zeta     New-to-used transition rate', zb, zc, '%11.3f', '%15.3f', ...
    sprintf('New-to-used capital ratio: %.3f', TGT.KnKu));

% kappa exists only where used capital does.
kb = getf(E(ib).fit,'kappa');
kc = NaN; if ~isempty(ic) && has_used(E(ic)), kc = getf(E(ic).fit,'kappa'); end
prow3('kappa    Breakage rate of used capital', kb, kc, '%11.3f', '%15.3f', ...
    sprintf('Used-share gradient: %.3f', TGT.su));

fprintf('\n  Notes: All targets are computed on the Vietnamese Small and Medium\n');
fprintf('  Enterprise Survey (CIEM et al., 2005-2015) and are listed in\n');
fprintf('  paper_targets.m. The two economies are calibrated separately. A dash\n');
fprintf('  marks a parameter that is undefined in an economy with no second-hand\n');
fprintf('  market, not one that is zero.\n');


%% ---------------------------------------------------------------- TABLE 4
head(L, 'TABLE 4.  Used Capital and Capital Misallocation');
if isempty(ic)
    fprintf('\n  (the model without used capital is missing; Table 4 needs it)\n');
else
    fprintf('\n  %-42s%12s%16s%25s\n', '', 'Benchmark', 'No used capital', 'Effect of used capital');
    fprintf('\n');
    erow('TFP loss, M (log points)',   E(ib).loss,               E(ic).loss,               '%12.2f', 'elim',  'eliminated');
    erow('Mean scale relative to target, K/K*',  getf(E(ib).base,'kk_mean'), getf(E(ic).base,'kk_mean'), '%12.3f', 'ratio', 'higher');
    erow('sd(log MRPK)',                         getana(E(ib),'sd_logMRPK'), getana(E(ic),'sd_logMRPK'), '%12.3f', 'elim',  'lower');
    erow('Share of firms constrained',           getf(E(ib).base,'constr'),  getf(E(ic).base,'constr'),  '%12.3f', 'elim',  'lower');

    fprintf('\n  Notes: The benchmark contains new and used capital, while the other\n');
    fprintf('  model is without used capital and recalibrated to the common targets.\n');
    fprintf('  Each economy is compared with its own frictionless twin (theta -> 1e5),\n');
    fprintf('  so the misallocation loss is a within-economy statistic:\n');
    fprintf('  M = 100 log(TFP frictionless / TFP constrained).\n');
end


%% ---------------------------------------------------------------- TABLE 5
head(L, 'TABLE 5.  Used Capital and Entry by Productive Poor Agents');
okS = ~isempty(ic) && has_sel(E(ib)) && has_sel(E(ic));
if ~okS
    fprintf('\n  (needs the life-cycle analyses; run analysis_used_capital.m and\n');
    fprintf('   analysis_cf.m in their folders first)\n');
else
    fprintf('\n  %-42s%12s%16s%25s\n', '', 'Benchmark', 'No used capital', 'Effect of used capital');
    fprintf('\n');
    pb = E(ib).ana.econ(1).pt_rate;   pc = E(ic).ana.econ(1).pt_rate;
    rb = E(ib).ana.econ(1).ptr_rich;  rc = E(ic).ana.econ(1).ptr_rich;
    fb = E(ib).ana.econ(2).pt_rate;   fc = E(ic).ana.econ(2).pt_rate;
    qb = E(ib).ana.econ(2).ptr_rich;  qc = E(ic).ana.econ(2).ptr_rich;

    pprow('Poor talented entrepreneurship rate', pb, pc);
    pprow('Rich talented entrepreneurship rate', rb, rc);
    fprintf('\n');
    dashrow('Poor talented rate, frictionless',  fb, fc);
    dashrow('Rich talented rate, frictionless',  qb, qc);
    fprintf('\n');

    Gb = 100*log(fb/pb);   Gc = 100*log(fc/pc);
    fprintf('  %-42s%12.1f%16.1f', 'Entry shortfall, S (log points)', Gb, Gc);
    if isfinite(Gb) && isfinite(Gc) && abs(Gc) > 1e-12
        fprintf('%22.1f%% repaired\n', 100*(1 - Gb/Gc));
    else
        fprintf('%25s\n', '--');
    end

    fprintf('\n  Notes: "Talented" denotes the top productivity decile and "poor"\n');
    fprintf('  denotes wealth below the benchmark population median. These cutoffs\n');
    fprintf('  are constructed once in the benchmark and applied in levels to both\n');
    fprintf('  economies. Each economy is compared with its own frictionless\n');
    fprintf('  counterpart, so the frictionless rows differ across columns. The\n');
    fprintf('  entry shortfall S is defined for poor talented agents.\n');
end


%% ---------------------------------------------------------------- TABLE 6
head(L, 'TABLE 6.  The result under three designs without used capital');
fprintf('\n  %-32s%8s%10s%15s%11s%9s%13s\n', 'Economy', 'eta_p', 'TFP loss', 'Misallocation', 'Entry', 'K/K*', 'sd(logMRPK)');
fprintf('  %-32s%8s%10s%15s%11s\n', '', '', '(log pts)', 'reduced', 'repaired');
fprintf('\n');

Gb = NaN;
if has_sel(E(ib)), Gb = 100*log(E(ib).ana.econ(2).pt_rate / E(ib).ana.econ(1).pt_rate); end

ORDER = { 'Benchmark (used K)',            ib
          'No used capital, tail pinned',  ip
          'Tail diagnostic',               id
          'No used capital, recalibrated', ic };
for j = 1:size(ORDER,1)
    i = ORDER{j,2};
    if isempty(i), continue; end
    fprintf('  %-32s', ORDER{j,1});
    fnum(getf(E(i).fit,'eta_p'), '%8.3f');
    fnum(E(i).loss,              '%10.2f');
    if i == ib
        fprintf('%15s%11s', '--', '--');
    else
        pct(100*(1 - E(ib).loss/E(i).loss), '%14.1f%%');
        en = NaN;
        if has_sel(E(i)) && isfinite(Gb)
            en = 100*(1 - Gb/(100*log(E(i).ana.econ(2).pt_rate / E(i).ana.econ(1).pt_rate)));
        end
        pct(en, '%10.1f%%');
    end
    fnum(getf(E(i).base,'kk_mean'),  '%9.3f');
    fnum(getana(E(i),'sd_logMRPK'),  '%13.3f');
    fprintf('\n');
end
fprintf('\n  Notes: The ladder is ordered by how much is allowed to re-adjust.\n');
fprintf('  "Tail pinned" holds the ability process at its benchmark value and\n');
fprintf('  recalibrates only theta and beta. "Tail diagnostic" begins from the\n');
fprintf('  recalibrated economy and returns eta_p to its benchmark value without\n');
fprintf('  refitting any moment; its own moments therefore drift, and it is a\n');
fprintf('  diagnostic rather than a calibrated column. Misallocation reduced is\n');
fprintf('  1 - M_benchmark/M_j; entry repaired is 1 - S_benchmark/S_j.\n');


%% ---------------------------------------------------------------- TABLE 7
head(L, 'TABLE 7.  Sensitivity to the two assigned CES parameters');
RPF = fullfile(HERE, 'benchmark', 'robustness_phi_result.mat');
RGF = fullfile(HERE, 'benchmark', 'robustness_gamma_result.mat');
haveP = exist(RPF, 'file') == 2;
haveG = exist(RGF, 'file') == 2;
if ~haveP && ~haveG
    fprintf('\n  (no sensitivity results found; run benchmark/robustness_phi.m\n');
    fprintf('   and benchmark/robustness_gamma.m)\n');
else
    P = []; G = []; g1 = []; g2 = [];
    if haveP
        Rp = load(RPF); rn = fieldnames(Rp); P = Rp.(rn{1});
        SRC(end+1).file = 'benchmark/robustness_phi_result.mat'; %#ok<AGROW>
        SRC(end).what   = 'Table 7, phi columns';
    end
    if haveG
        Rg = load(RGF); rn = fieldnames(Rg); G = Rg.(rn{1});
        Sg = G.S;
        if numel(Sg) >= 2, g1 = Sg(1); g2 = Sg(2); end
        SRC(end+1).file = 'benchmark/robustness_gamma_result.mat'; %#ok<AGROW>
        SRC(end).what   = 'Table 7, gamma columns';
    end
    if ~haveP, fprintf('\n  (phi columns absent; run benchmark/robustness_phi.m)\n');  end
    if ~haveG, fprintf('\n  (gamma columns absent; run benchmark/robustness_gamma.m)\n'); end

    % Two sweeps, one table. They are printed side by side but they are NOT a
    % single four-point grid: each sweep re-clears prices around the benchmark
    % calibration holding everything else fixed, and the two were solved in
    % separate runs. Read each sweep against its own reference, which is why
    % the phi block carries the benchmark column and the gamma block does not.
    fprintf('\n  %-44s%22s%22s\n', '', 'phi sweep', 'gamma sweep');
    fprintf('  %-44s%11s%11s%11s%11s\n', '', 'phi=0.70', 'phi=0.75', 'gamma=1.7', 'gamma=2.5');
    fprintf('  %-44s%11s\n', '', '(benchmark)');

    fprintf('\n  Panel A: Gain from removing the collateral constraint\n');
    qrow('TFP (%)', pf(P,'gains070','TFP'), pf(P,'gains075','TFP'), ...
                    gtfp(g1), gtfp(g2), '%11.2f');
    qrow('TFP loss M (log points)', pf(P,'gains070','L'), pf(P,'gains075','L'), ...
                    gf(g1,'gain_L'), gf(g2,'gain_L'), '%11.2f');
    % Averted is measured against the same no-used-capital economy in every
    % column, which is what makes the four comparable at all.
    Lc = NaN; if ~isempty(ic), Lc = E(ic).loss; end
    qrow('Misallocation reduced (%)', avert(pf(P,'gains070','L'), Lc), ...
                    avert(pf(P,'gains075','L'), Lc), ...
                    avert(gf(g1,'gain_L'), Lc), avert(gf(g2,'gain_L'), Lc), '%11.1f');
    qrow('Aggregate output (%)', NaN, NaN, gf(g1,'gain_Y'), gf(g2,'gain_Y'), '%11.2f');

    fprintf('\n  Panel B: Allocation in the constrained equilibrium\n');
    qrow('Mean scale relative to target, K/K*', pf(P,'f070','kk_mean'), pf(P,'f075','kk_mean'), ...
                    gfr(g1,'kk_mean'), gfr(g2,'kk_mean'), '%11.3f');
    qrow('Share of firms constrained', pf(P,'f070','constr'), pf(P,'f075','constr'), ...
                    gfr(g1,'constr'), gfr(g2,'constr'), '%11.3f');
    qrow('Entrepreneur share', pf(P,'f070','ent'), pf(P,'f075','ent'), ...
                    gfr(g1,'ent'), gfr(g2,'ent'), '%11.3f');

    fprintf('\n  Panel C: Cleared prices\n');
    qrow('Wage W',                pf(P,'f070','W'),    pf(P,'f075','W'),    gfr(g1,'W'),    gfr(g2,'W'),    '%11.3f');
    qrow('Interest rate r',       pf(P,'f070','r'),    pf(P,'f075','r'),    gfr(g1,'r'),    gfr(g2,'r'),    '%11.3f');
    qrow('New-capital rental Rn', pf(P,'f070','Rnew'), pf(P,'f075','Rnew'), gfr(g1,'Rnew'), gfr(g2,'Rnew'), '%11.3f');
    qrow('Used-capital price pu', pf(P,'f070','pu'),   pf(P,'f075','pu'),   gfr(g1,'pu'),   gfr(g2,'pu'),   '%11.3f');

    fprintf('\n  Panel D: Targeted moments\n');
    qrow('Debt-to-output ratio',        pm(P,'f070','DtoY'),      pm(P,'f075','DtoY'),      gm(g1,'DtoY'),      gm(g2,'DtoY'),      '%11.3f');
    qrow('Top-decile employment share', pm(P,'f070','top10'),     pm(P,'f075','top10'),     gm(g1,'top10'),     gm(g2,'top10'),     '%11.3f');
    qrow('Firm exit rate',              pm(P,'f070','exit'),      pm(P,'f075','exit'),      gm(g1,'exit'),      gm(g2,'exit'),      '%11.3f');
    qrow('Wage Gini',                   pm(P,'f070','wGini'),     pm(P,'f075','wGini'),     gm(g1,'wGini'),     gm(g2,'wGini'),     '%11.3f');
    qrow('Used-share gradient',         pm(P,'f070','beta_size'), pm(P,'f075','beta_size'), gm(g1,'beta_size'), gm(g2,'beta_size'), '%11.3f');
    qrow('New-to-used capital ratio',   pm(P,'f070','KnKu'),      pm(P,'f075','KnKu'),      gm(g1,'KnKu'),      gm(g2,'KnKu'),      '%11.3f');

    fprintf('\n  Notes: phi and gamma are assigned rather than calibrated, so each\n');
    fprintf('  sweep re-clears prices at the new value with every other parameter\n');
    fprintf('  held at its calibrated level. Panel D shows that the targeted\n');
    fprintf('  moments barely move, so the headline in Panel A does not rest on\n');
    fprintf('  either choice.\n');
    fprintf('\n  Read each sweep against itself. The two blocks were solved in\n');
    fprintf('  separate runs, so a level difference common to both gamma columns\n');
    fprintf('  is a run effect and not a gamma effect; the gamma result is the\n');
    fprintf('  difference BETWEEN 1.7 and 2.5. Across that range the misallocation\n');
    fprintf('  loss moves by about 0.06 log points, roughly one percent of its\n');
    fprintf('  level. Full discussion in benchmark/robustness_gamma_results.txt.\n');
end


%% -------------------------------------------------------------- PROVENANCE
head(L, 'APPENDIX.  Provenance of every number above');
fprintf('\n  Tables 1 and 2 are echoed from a committed results file. Tables 3 to 7\n');
fprintf('  are computed in tables.m from the files below. Nothing is hardcoded.\n');
fprintf('\n  %-52s%12s  %s\n', 'source file', 'bytes', 'last modified');
fprintf('  %s\n', repmat('-', 1, 76));
for i = 1:numel(SRC)
    d = dir(fullfile(HERE, SRC(i).file));
    if isempty(d)
        fprintf('  %-52s%12s  %s\n', SRC(i).file, 'MISSING', '--');
    else
        fprintf('  %-52s%12d  %s\n', SRC(i).file, d.bytes, datestr(d.datenum, 'yyyy-mm-dd HH:MM'));
    end
end

fprintf('\n  Field behind each row\n');
fprintf('  %s\n', repmat('-', 1, 76));
fprintf('    Table 3 Panel A   parameters_*.m, read as text (never executed)\n');
fprintf('    Table 3 Panel B   out.theta, out.eta_p, out.omega, out.sigma_y,\n');
fprintf('                      out.beta, out.kappa in *_result.mat;\n');
fprintf('                      zeta = (delta + kappa) / (Kn/Ku)\n');
fprintf('    Table 3 targets   paper_targets.m\n');
fprintf('    Table 4 M         100 log(nofriction_result.frictionless.TFP\n');
fprintf('                                / nofriction_result.base.TFP)\n');
fprintf('    Table 4 K/K*      nofriction_result.base.kk_mean\n');
fprintf('    Table 4 sd MRPK   analysis_*_result.A.econ(1).sd_logMRPK\n');
fprintf('    Table 4 constr    nofriction_result.base.constr\n');
fprintf('    Table 5 rates     analysis_*_result.A.econ(1).pt_rate, .ptr_rich\n');
fprintf('                      frictionless rows use A.econ(2)\n');
fprintf('    Table 5 S         100 log(econ(2).pt_rate / econ(1).pt_rate)\n');
fprintf('    Table 6           the same fields, one row per economy\n');
fprintf('    Table 7 phi       robustness_phi_result.{gains070,gains075,f070,f075}\n');
fprintf('    Table 7 gamma     robustness_gamma_result.S(i).{fric,nofr,gain_L,gain_Y}\n');

fprintf('\n  To re-derive the stored equilibria rather than read them, run\n');
fprintf('  benchmark/verify_benchmark.m. It ignores every cache and re-solves\n');
fprintf('  the benchmark from the parameters in parameters_benchmark.m, which\n');
fprintf('  takes about forty minutes on the production grid.\n');

L('=');
fprintf('  Written to %s\n', OUT);
L('=');
diary off;
end


%% ======================================================================
%  helpers
%  ======================================================================

function [E, SRC] = load_all(SPEC)
E = struct('name',{},'fit',{},'base',{},'fri',{},'loss',{},'ana',{});
SRC = struct('file',{},'what',{});
for i = 1:size(SPEC,1)
    d   = SPEC{i,2};
    fnf = fullfile(d, 'nofriction_result.mat');
    ffi = fullfile(d, SPEC{i,3});
    if ~exist(fnf,'file') || ~exist(ffi,'file'), continue; end
    Rn = load(fnf); rn = fieldnames(Rn); Rn = Rn.(rn{1});
    Qf = load(ffi);
    j = numel(E) + 1;
    E(j).name = SPEC{i,1};
    E(j).fit  = Qf.out;
    E(j).base = Rn.base;
    E(j).fri  = Rn.frictionless;
    E(j).loss = 100 * log(Rn.frictionless.TFP / Rn.base.TFP);
    SRC(end+1).file = strrep(ffi, filesep, '/'); %#ok<AGROW>
    SRC(end).what   = SPEC{i,1};
    SRC(end+1).file = strrep(fnf, filesep, '/'); %#ok<AGROW>
    SRC(end).what   = [SPEC{i,1} ' frictionless twin'];
    fana = fullfile(d, SPEC{i,4});
    if exist(fana,'file')
        Aa = load(fana); E(j).ana = Aa.A;
        SRC(end+1).file = strrep(fana, filesep, '/'); %#ok<AGROW>
        SRC(end).what   = [SPEC{i,1} ' selection and dispersion'];
    else
        E(j).ana = [];
    end
    % A fit file and a solve that are not the same economy would misreport
    % every column, so say so loudly rather than print it.
    if isfield(E(j).fit,'theta') && isfield(E(j).base,'theta') && ...
            abs(E(j).fit.theta - E(j).base.theta) > 1e-6
        warning('tables:mismatch', ...
            '%s: fit theta %.6f but solved theta %.6f -- not the same economy.', ...
            SPEC{i,1}, E(j).fit.theta, E(j).base.theta);
    end
end
end

function P = read_params(pfile)
% Read assigned parameters out of a parameters_*.m file. The file is parsed as
% text and never executed, so no global is touched and nothing can be
% overwritten by running this.
P = struct();
if ~exist(pfile,'file'), return; end
txt = fileread(pfile);
NAMES = {'alpha','delta','eta_span','sigma_crra','rho_y','phi','gamma'};
for j = 1:numel(NAMES)
    tok = regexp(txt, ['^\s*' NAMES{j} '\s*=\s*([^;%]+);'], 'tokens', 'lineanchors', 'once');
    if isempty(tok), continue; end
    s = strtrim(tok{1});
    v = str2double(s);
    if ~isfinite(v)
        % Assigned constants are sometimes written as a fraction, e.g. 2/3.
        v = local_fraction(s);
    end
    if isfinite(v), P.(NAMES{j}) = v; end
end
end

function v = local_fraction(s)
v = NaN;
t = regexp(s, '^\s*([\d.]+)\s*/\s*([\d.]+)\s*$', 'tokens', 'once');
if ~isempty(t)
    a = str2double(t{1}); b = str2double(t{2});
    if isfinite(a) && isfinite(b) && b ~= 0, v = a/b; end
end
end

function v = gp(P, f)
if isstruct(P) && isfield(P, f), v = P.(f); else, v = NaN; end
end

function v = getf(s, f)
if isstruct(s) && isfield(s, f) && isfinite(s.(f)), v = s.(f); else, v = NaN; end
end

function v = getana(Ei, f)
v = NaN;
if ~isempty(Ei.ana) && isfield(Ei.ana.econ(1), f) && isfinite(Ei.ana.econ(1).(f))
    v = Ei.ana.econ(1).(f);
end
end

function tf = has_sel(Ei)
tf = ~isempty(Ei.ana) && isfield(Ei.ana.econ(1), 'pt_rate') && isfield(Ei.ana.econ(1), 'ptr_rich');
end

function tf = has_used(Ei)
% An economy has a second-hand market only if used capital is actually priced.
k = getf(Ei.fit, 'kappa');
tf = isfinite(k) && k > 0;
end

function z = zeta_of(E, i, TGT)
% zeta = (delta + kappa) / (Kn/Ku). With no used capital there is no
% transition at all, so zeta is zero rather than undefined.
z = NaN;
if isempty(i), return; end
if ~has_used(E(i)), z = 0; return; end
kap = getf(E(i).fit, 'kappa');
if isfinite(kap) && isfinite(TGT.KnKu) && TGT.KnKu > 0
    z = (0.06 + kap) / TGT.KnKu;
end
end

function head(L, txt)
fprintf('\n\n'); L('='); fprintf('  %s\n', txt); L('=');
end

function fnum(v, fmt)
if isfinite(v), fprintf(fmt, v); else, fprintf(strrep_width(fmt), '--'); end
end

function f = strrep_width(fmt)
w = regexp(fmt, '%(\d+)', 'tokens', 'once');
if isempty(w), f = '%s'; else, f = ['%' w{1} 's']; end
end

function pct(v, fmt)
if isfinite(v), fprintf(fmt, v); else, fprintf(strrep_width(fmt), '--'); end
end

function arow(lab, vb, vc, src)
% One Panel A row of Table 3.
fprintf('  %-44s', lab);
fnum(vb, '%11.3f'); fnum(vc, '%15.3f');
fprintf('   %s\n', src);
end

function brow(E, ib, ic, fld, lab, src)
% One Panel B row of Table 3, read from the solved economies.
vb = NaN; vc = NaN;
if ~isempty(ib), vb = getf(E(ib).fit, fld); end
if ~isempty(ic), vc = getf(E(ic).fit, fld); end
prow3(lab, vb, vc, '%11.3f', '%15.3f', src);
end

function prow3(lab, vb, vc, fb, fc, src)
fprintf('  %-44s', lab);
fnum(vb, fb); fnum(vc, fc);
fprintf('   %s\n', src);
end

function erow(lab, vb, vc, fmt, mode, word)
% One Table 4 row. MODE 'elim' reports 1 - benchmark/other; 'ratio' reports
% benchmark/other - 1. WORD is the noun the paper attaches to it.
fprintf('  %-42s', lab);
fnum(vb, fmt); fnum(vc, fmt);
if isfinite(vb) && isfinite(vc) && abs(vc) > 1e-12
    switch mode
        case 'elim',  p = 100*(1 - vb/vc);
        case 'ratio', p = 100*(vb/vc - 1);
        otherwise,    p = NaN;
    end
    w = word;
    if strcmp(mode,'ratio') && p < 0, w = 'lower'; p = -p; end
    fprintf('%16.1f%% %s\n', p, w);
else
    fprintf('%25s\n', '--');
end
end

function pprow(lab, vb, vc)
% Table 5 entry-rate row: the effect is a difference in percentage points.
fprintf('  %-42s', lab);
fnum(vb, '%12.3f'); fnum(vc, '%16.3f');
if isfinite(vb) && isfinite(vc)
    fprintf('%+20.1f p.p.\n', 100*(vb - vc));
else
    fprintf('%25s\n', '--');
end
end

function dashrow(lab, vb, vc)
fprintf('  %-42s', lab);
fnum(vb, '%12.3f'); fnum(vc, '%16.3f');
fprintf('%25s\n', '--');
end

function trow(lab, v70, v75, fmt)
fprintf('  %-44s', lab);
fnum(v70, fmt); fnum(v75, fmt);
fprintf('\n');
end

function echo_empirical(L, tfile)
% Tables 1 and 2 come from the survey microdata, which is not included here.
% They are echoed from the committed results file the Python script wrote, so
% what appears here is exactly what that script produced.
head(L, 'TABLES 1 AND 2.  Used capital in Vietnamese SMEs (empirical)');
if ~exist(tfile, 'file')
    fprintf('\n  %s is missing.\n', tfile);
    fprintf('  Produce it with:  cd empirical/motivation && python used_capital_table.py\n');
    fprintf('  That needs the survey rounds in empirical/enterprise_data/, which are\n');
    fprintf('  not included here. Download them from\n');
    fprintf('    https://www.wider.unu.edu/database/viet-nam-sme-database\n');
    fprintf('  See empirical/README.md.\n');
    return
end
fprintf('\n  Echoed verbatim from empirical/motivation/used_capital_table.txt,\n');
fprintf('  written by used_capital_table.py from the SME survey. That microdata\n');
fprintf('  is not included in this package, so the committed results file is\n');
fprintf('  the audit trail; the full notes are in it.\n\n');
fid = fopen(tfile, 'r');
if fid < 0, fprintf('  (could not open %s)\n', tfile); return; end
started = false;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if ~started && contains(ln, 'TABLE 1.'), started = true; end
    if started
        if ~isempty(strtrim(ln)) && contains(ln, 'Notes') && ~contains(ln, 'TABLE')
            break
        end
        fprintf('%s\n', ln);
    end
end
fclose(fid);
end

function v = pf(P, blk, f)
% One field from the phi sweep, tolerant of an absent result file.
v = NaN;
if isempty(P) || ~isstruct(P) || ~isfield(P, blk), return; end
b = P.(blk);
if isstruct(b) && isfield(b, f) && isfinite(b.(f)), v = b.(f); end
end

function v = pm(P, blk, f)
% One targeted moment from the phi sweep.
v = NaN;
if isempty(P) || ~isstruct(P) || ~isfield(P, blk), return; end
b = P.(blk);
if ~isstruct(b) || ~isfield(b, 'm'), return; end
m = b.m;
if isstruct(m) && isfield(m, f) && isfinite(m.(f)), v = m.(f); end
end

function v = gf(Sg, f)
% One top-level field from one gamma point.
v = NaN;
if isempty(Sg) || ~isstruct(Sg) || ~isfield(Sg, f), return; end
if isfinite(Sg.(f)), v = Sg.(f); end
end

function v = gfr(Sg, f)
% One field of the constrained economy at one gamma point.
v = NaN;
if isempty(Sg) || ~isstruct(Sg) || ~isfield(Sg, 'fric'), return; end
fr = Sg.fric;
if isstruct(fr) && isfield(fr, f) && isfinite(fr.(f)), v = fr.(f); end
end

function v = gm(Sg, f)
% One targeted moment at one gamma point.
v = NaN;
if isempty(Sg) || ~isstruct(Sg) || ~isfield(Sg, 'fric'), return; end
fr = Sg.fric;
if ~isstruct(fr) || ~isfield(fr, 'm'), return; end
m = fr.m;
if isstruct(m) && isfield(m, f) && isfinite(m.(f)), v = m.(f); end
end

function v = gtfp(Sg)
% The gamma file stores the two economies rather than the TFP gain, so the
% gain is formed here the same way the phi file forms its own: the percentage
% difference between the frictionless twin and the constrained economy.
v = NaN;
if isempty(Sg) || ~isstruct(Sg), return; end
if ~isfield(Sg,'fric') || ~isfield(Sg,'nofr'), return; end
a = Sg.fric; b = Sg.nofr;
if isstruct(a) && isstruct(b) && isfield(a,'TFP') && isfield(b,'TFP') ...
        && isfinite(a.TFP) && isfinite(b.TFP) && a.TFP > 0
    v = 100*(b.TFP/a.TFP - 1);
end
end

function v = avert(Lj, Lc)
% Share of the no-used-capital misallocation loss that this economy avoids.
v = NaN;
if isfinite(Lj) && isfinite(Lc) && Lc > 0, v = 100*(1 - Lj/Lc); end
end

function qrow(lab, v1, v2, v3, v4, fmt)
fprintf('  %-44s', lab);
fnum(v1, fmt); fnum(v2, fmt); fnum(v3, fmt); fnum(v4, fmt);
fprintf('\n');
end
