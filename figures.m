function figures()
% figures -- the figures the paper reports, drawn in about a minute from
% results that are already stored in this package. Nothing is re-solved and
% no curve in this file is drawn by hand.
%
%   >> figures
%
% Writes fig_mechanism.pdf/.png, fig_mrpk_by_wealth.pdf/.png,
% fig_usedcapital_tilt.pdf/.png  (names match the manuscript)
% and a transcript in FIGURES.txt.
%
%   Figure 1  The mechanism of used capital                 two panels
%             left   shadow cost mu by benchmark entrepreneur wealth decile,
%                    for low, mid and high entrepreneurial productivity z
%             right  used-capital value share against mu, at benchmark prices
%   Figure 2  Used capital and the allocation of capital
%             MRPK relative to each economy's own frictionless level, by
%             benchmark entrepreneur wealth percentile
%
% WHERE THE CURVES COME FROM
%   Every point is read from the solved equilibria in
%   benchmark/analysis_used_capital_result.mat and
%   model_with_no_used/analysis_cf_result.mat. Those files hold the invariant
%   distribution over (a, z) together with each firm's shadow cost mu, its
%   used-capital share and its MRPK, so the curves are weighted averages over
%   the model's own population rather than anything typed in.
%
%   FIGURES.txt prints the numbers behind every curve -- the ten decile values
%   for each productivity group, the binned (mu, share) pairs, and the
%   percentile profile of MRPK in both economies -- so any point on any curve
%   can be checked against the picture by eye.
%
%   To re-derive the stored equilibria rather than read them, run
%   benchmark/verify_benchmark.m, then analysis_used_capital.m and
%   analysis_cf.m in their folders.

HERE = fileparts(mfilename('fullpath'));
cd(HERE); addpath(HERE);

OUT = fullfile(HERE, 'FIGURES.txt');
diary off;
try, if exist(OUT,'file'), delete(OUT); end, catch, end   %#ok<TRYNC>
diary(OUT); diary on;

L = @(c) fprintf('%s\n', repmat(c, 1, 78));
L('=');
fprintf('  FINANCIAL FRICTIONS, USED CAPITAL, AND MISALLOCATION\n');
fprintf('  Figures 1 and 2, drawn from stored results\n');
fprintf('  generated %s\n', char(datetime('now')));
L('=');

% ---- load the two economies ----------------------------------------------
BF = fullfile(HERE, 'benchmark',          'analysis_used_capital_result.mat');
CF = fullfile(HERE, 'model_with_no_used', 'analysis_cf_result.mat');
for f = {BF, CF}
    if ~exist(f{1}, 'file')
        diary off;
        error('figures:missingInput', ...
            ['%s is missing. Run analysis_used_capital.m in benchmark/ and ' ...
             'analysis_cf.m in model_with_no_used/ first.'], f{1});
    end
end
Bm = load(BF);  B = Bm.A;
Cm = load(CF);  C = Cm.A;
bf = B.econ(1);   bn = B.econ(2);      % benchmark: constrained, frictionless
cf = C.econ(1);   cn = C.econ(2);      % no used capital: same two

fprintf('\n  source files\n');
fprintf('  %s\n', repmat('-', 1, 76));
for f = {BF, CF}
    d = dir(f{1});
    fprintf('  %-52s%10d B  %s\n', ...
        strrep(erase(f{1}, [HERE filesep]), filesep, '/'), ...
        d.bytes, datestr(d.datenum, 'yyyy-mm-dd HH:MM'));
end

% ---- house style ---------------------------------------------------------
CB = [0.06 0.20 0.44];    % benchmark   (navy)
CC = [0.94 0.11 0.56];    % counterfactual (magenta)
CG = [0.45 0.45 0.45];    % efficient reference
LW = 1.7;
FS = 8;
set(groot, 'defaultAxesFontSize', FS);
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultAxesLineWidth', 0.75);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesTickLength', [0.012 0.012]);
set(groot, 'defaultAxesXColor', [0.25 0.25 0.25]);
set(groot, 'defaultAxesYColor', [0.25 0.25 0.25]);
set(groot, 'defaultAxesGridAlpha', 0.10);
set(groot, 'defaultAxesBox', 'off');


%% ===================================================================== FIG 1
% Two panels on one canvas: the shadow cost across the wealth distribution,
% then what firms do about it. Explicit axes positions rather than tiledlayout
% so the file runs on older MATLAB releases too.
f1  = newfig(7.4, 3.0);
axL = axes(f1, 'Position', [0.070 0.150 0.390 0.795]); hold(axL, 'on');
axR = axes(f1, 'Position', [0.580 0.150 0.390 0.795]); hold(axR, 'on');

% ---- left panel: mu by wealth decile, split by productivity ---------------
SH  = [0.45 0.24 0.00];                 % light -> dark
MKA = {'o', 's', '^'};
TLA = {'low \itz', 'mid \itz', 'high \itz'};
xdF = linspace(1, 10, 300);
hA  = gobjects(3,1);
for t = 1:3
    col = CB + (1 - CB)*SH(t);
    plot(axL, xdF, interp1((1:10)', bf.tabMU10(:,t), xdF, 'pchip'), '-', ...
        'Color', col, 'LineWidth', LW, 'HandleVisibility', 'off');
    hA(t) = plot(axL, (1:10)', bf.tabMU10(:,t), MKA{t}, 'Color', col, ...
        'LineStyle', 'none', 'MarkerFaceColor', col, 'MarkerSize', 4.2);
end
xlabel(axL, 'benchmark entrepreneur wealth decile');
ylabel(axL, 'shadow cost \mu');
xlim(axL, [0.5 10.5]); xticks(axL, 1:10);
yl = ylim(axL); ylim(axL, [0 yl(2)]);
legend(axL, hA, TLA, 'Location', 'northeast', 'Box', 'off', 'FontSize', FS-1);
grid(axL, 'on');

fprintf('\n\n  FIGURE 1, LEFT PANEL -- shadow cost mu by wealth decile\n');
fprintf('  source: A.econ(1).tabMU10, a 10-by-3 table in the benchmark file\n\n');
fprintf('    %-8s%12s%12s%12s\n', 'decile', 'low z', 'mid z', 'high z');
for i = 1:10
    fprintf('    %-8d%12.4f%12.4f%12.4f\n', i, bf.tabMU10(i,1), bf.tabMU10(i,2), bf.tabMU10(i,3));
end

% ---- right panel: the used-capital value share against mu ----------------
% The share is a value share, p_u K_u / (p_n K_n + p_u K_u), with p_n = 1.
% Firms are binned by mu over the invariant distribution and averaged within
% bin, so each plotted point is a population average, not a fitted curve.
mu_b = bf.mu(:);  s_b = bf.ushare(:);  w_b = bf.me(:)/sum(bf.me);
isU  = mu_b <= 1e-8;                 % the unconstrained atom, plotted apart
isC  = ~isU;
muC  = mu_b(isC);  sC = s_b(isC);  wC = w_b(isC);

NBIN = 20;
edm  = arrayfun(@(q) wquantile_vec(muC, wC, q), linspace(0, 1, NBIN+1));
edm(1) = -inf;  edm(end) = inf;
xm = nan(1,NBIN);  ys = nan(1,NBIN);
for i = 1:NBIN
    k = muC >= edm(i) & muC < edm(i+1);
    if ~any(k), continue; end
    ww = wC(k)/sum(wC(k));
    xm(i) = sum(muC(k).*ww);  ys(i) = sum(sC(k).*ww);
end
ok = ~isnan(xm);  xm = xm(ok);  ys = ys(ok);

s_unc = sum(s_b(isU).*w_b(isU)) / sum(w_b(isU));

% The top mu bin sits an order of magnitude above the bottom one and would
% compress everything else into the left edge, so the axis stops at the 95th
% percentile of mu among constrained firms. The dropped bins are listed below.
XCUT = wquantile_vec(muC, wC, 0.95);
keep = xm <= XCUT;

% The same panel is drawn twice: once as the right half of Figure 1 and once on
% its own canvas, so the mechanism figure can be used either way. Both come from
% one computation, so the two cannot disagree.
y_unc = 100*s_unc;
y_end = 100*ys(find(keep, 1, 'last'));

draw_tilt(axR, xm, ys, keep, s_unc, XCUT, y_unc, y_end, CB, LW, FS);

f1b  = onepanel();
ax1b = axes(f1b); hold(ax1b, 'on');
draw_tilt(ax1b, xm, ys, keep, s_unc, XCUT, y_unc, y_end, CB, LW, FS);

fprintf('\n\n  FIGURE 1, RIGHT PANEL -- used-capital value share against mu\n');
fprintf('  source: A.econ(1).mu, .ushare and .me; firms binned by mu into %d\n', NBIN);
fprintf('  equal-measure bins, averaged within bin using the invariant measure\n\n');
fprintf('    %-14s%14s%14s\n', 'bin', 'mean mu', 'share (%)');
fprintf('    %-14s%14.4f%14.2f\n', 'unconstrained', 0, y_unc);
for i = 1:numel(xm)
    tag = sprintf('%d', i);
    if ~keep(i), tag = sprintf('%d (off axis)', i); end
    fprintf('    %-14s%14.4f%14.2f\n', tag, xm(i), 100*ys(i));
end
fprintf('\n    unconstrained share      %6.2f %%\n', y_unc);
fprintf('    at the axis cut          %6.2f %%   (mu = %.4f)\n', y_end, XCUT);
fprintf('    rise plotted             %+6.2f pp\n', y_end - y_unc);
fprintf('    mass beyond the cut      %6.2f %% of constrained firms\n', ...
    100*sum(wC(muC > XCUT))/sum(wC));

savefig_paper(f1,  'fig_mechanism', 'both');
savefig_paper(f1b, 'fig_usedcapital_tilt', 'both');
fprintf('\n  [saved] fig_mechanism.pdf and fig_mechanism.png\n');
fprintf('  [saved] fig_usedcapital_tilt.pdf and .png  (right panel, standalone)\n');


%% ===================================================================== FIG 2
% MRPK relative to each economy's OWN frictionless level. Dividing each
% economy by its own efficient level is what makes the two comparable: the
% no-used-capital economy is recalibrated, so its levels differ.
mrpk_eff   = wmean(bn.MRPK, bn.me);
mrpk_eff_c = wmean(cn.MRPK, cn.me);

pct    = linspace(2, 98, 21);
edges1 = arrayfun(@(p) weighted_quantile(bf.a, bf.me, p/100), pct);
mass_c = bin_mass(cf.a, cf.me, edges1);
trimmed = false;
if min(mass_c) < 0.005 * sum(cf.me)
    % A bin holding under half a percent of the counterfactual's entrepreneurs
    % would plot a mean over almost nobody, so the range is trimmed instead.
    pct     = linspace(5, 98, 20);
    trimmed = true;
end

% Ability-standardised: each wealth bin is averaged over a COMMON productivity
% distribution, so the profile is not contaminated by the two economies
% sorting productivity differently across wealth.
[wb, mb] = profile_by_wealth_std(bf.a, bf.MRPK, bf.z, bf.me, bf.a, bf.me, pct, 'median');
[~,  mc] = profile_by_wealth_std(cf.a, cf.MRPK, cf.z, cf.me, bf.a, bf.me, pct, 'median');
mbn = mb / mrpk_eff;   mcn = mc / mrpk_eff_c;

f2 = onepanel();  ax = axes(f2); hold(ax, 'on');
xf = linspace(min(wb), max(wb), 400);
h_cf = plot(ax, xf, interp1(wb, mcn, xf, 'pchip'), '--', 'Color', CC, 'LineWidth', LW);
h_bm = plot(ax, xf, interp1(wb, mbn, xf, 'pchip'), '-',  'Color', CB, 'LineWidth', LW);
yline(ax, 1, '--', 'Color', CG, 'LineWidth', 1.1);
text(ax, 4, 1, 'efficient ', 'Color', CG, 'FontSize', FS-1, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
xlabel(ax, 'Benchmark entrepreneur wealth percentile');
ylabel(ax, 'MRPK / MRPK^*');
xlim(ax, [0 100]); xticks(ax, 0:20:100);
tighten_y(ax, 'floor', 0.9);
legend(ax, [h_bm h_cf], {'Benchmark (with used capital)', 'Without used capital'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', FS-1, ...
    'TextColor', [0.15 0.15 0.15]);
grid(ax, 'on');

fprintf('\n\n  FIGURE 2 -- MRPK relative to each economy''s own frictionless level\n');
fprintf('  source: A.econ(1).MRPK, .a, .z, .me in both files;\n');
fprintf('  the efficient level is the measure-weighted mean of A.econ(2).MRPK\n\n');
fprintf('    efficient MRPK, benchmark          %10.4f\n', mrpk_eff);
fprintf('    efficient MRPK, no used capital    %10.4f\n', mrpk_eff_c);
if trimmed
    fprintf('    percentile range trimmed to [5, 98]: a bin held under 0.5%% of the\n');
    fprintf('    no-used-capital entrepreneurs\n');
end
fprintf('\n    %-14s%16s%18s\n', 'percentile', 'benchmark', 'no used capital');
for i = 1:numel(wb)
    fprintf('    %-14.1f%16.3f%18.3f\n', wb(i), mbn(i), mcn(i));
end

savefig_paper(f2, 'fig_mrpk_by_wealth', 'both');
fprintf('\n  [saved] fig_mrpk_by_wealth.pdf and .png\n');

L('=');
fprintf('  Transcript written to %s\n', OUT);
fprintf('  Every plotted point is listed above and can be checked against the\n');
fprintf('  picture by eye. Nothing in this file is a hardcoded coordinate.\n');
L('=');
diary off;
end


%% ======================================================================
%  helpers -- copied unchanged from make_figures.m so the output matches
%  ======================================================================

function f = newfig(win, hin)
% Journal-width canvas in inches, so fonts land at the right size after LaTeX
% includes the file at \textwidth.
f = figure('Visible','off', 'Units','inches', ...
    'Position', [1 1 win hin], 'Color', 'w', ...
    'PaperUnits','inches', 'PaperSize', [win hin], ...
    'PaperPositionMode','auto');
end

function f = onepanel()
f = newfig(3.6, 3.0);
end

function savefig_paper(f, stem, fmt)
switch lower(fmt)
    case 'pdf'
        exportgraphics(f, [stem '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    case 'png'
        exportgraphics(f, [stem '.png'], 'Resolution', 300, 'BackgroundColor','white');
    case 'both'
        exportgraphics(f, [stem '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        exportgraphics(f, [stem '.png'], 'Resolution', 300, 'BackgroundColor','white');
    otherwise
        error('savefig_paper: unknown format ''%s''.', fmt);
end
close(f);
end

function tighten_y(ax, varargin)
% Shrink the y-axis to the data plus a small pad instead of leaving the large
% empty band MATLAB's autoscale often produces.
p = inputParser;
addParameter(p, 'pad', 0.06);
addParameter(p, 'floor', -inf);
addParameter(p, 'ceil', inf);
parse(p, varargin{:});
ch = findobj(ax, 'Type', 'line');
v = [];
for k = 1:numel(ch)
    d = get(ch(k), 'YData');
    v = [v, d(isfinite(d))]; %#ok<AGROW>
end
if isempty(v), return; end
lo = min(v); hi = max(v); rg = max(hi - lo, eps);
lo = max(lo - p.Results.pad*rg, p.Results.floor);
hi = min(hi + p.Results.pad*rg, p.Results.ceil);
if hi > lo, ylim(ax, [lo hi]); end
end

function m = bin_mass(a, w, edges)
a = a(:); w = w(:);
ok = isfinite(a) & isfinite(w) & w > 0;
m  = zeros(1, numel(edges)-1);
for i = 1:numel(edges)-1
    m(i) = sum(w(ok & a >= edges(i) & a < edges(i+1)));
end
end

function q = wquantile_vec(v, w, p)
% Weighted quantiles at every p in one sort, rather than one sort per p.
v = v(:); w = w(:);
ok = isfinite(v) & isfinite(w) & w > 0;
[vs, ix] = sort(v(ok));  ws = w(ok);  ws = ws(ix);
c = cumsum(ws) / sum(ws);
[cu, iu] = unique(c, 'last');
q = interp1(cu, vs(iu), p(:).', 'linear', 'extrap');
end

function q = weighted_quantile(v, w, p)
v = v(:); w = w(:);
ok = isfinite(v) & isfinite(w) & w > 0;
[vs, ix] = sort(v(ok));  ws = w(ok);  ws = ws(ix);
c = cumsum(ws) / sum(ws);
k = find(c >= p, 1, 'first');
if isempty(k), k = numel(vs); end
q = vs(k);
end

function [xc, yc] = profile_by_wealth_std(a, y, z, w, aref, wref, pct, stat)
% Statistic of y in common wealth-percentile bins, directly standardised on
% the productivity distribution so the two economies are comparable.
if nargin < 8 || isempty(stat), stat = 'mean'; end
a = a(:); y = y(:); z = z(:); w = w(:);
ok = isfinite(a) & isfinite(y) & isfinite(w) & w > 0;
edges = arrayfun(@(p) weighted_quantile(aref, wref, p/100), pct);
zs = unique(z(ok));
pz = arrayfun(@(v) sum(w(ok & z == v)), zs);
pz = pz / max(sum(pz), 1e-300);
xc = nan(1, numel(edges)-1);  yc = nan(1, numel(edges)-1);
for i = 1:numel(edges)-1
    if edges(i+1) <= edges(i), continue; end
    sel = ok & a >= edges(i) & a < edges(i+1);
    xc(i) = 0.5*(pct(i) + pct(i+1));
    if ~any(sel), continue; end
    num = 0; den = 0;
    for j = 1:numel(zs)
        cell_j = sel & z == zs(j);
        wj = sum(w(cell_j));
        if wj > 0
            if strcmpi(stat, 'median')
                cell_val = weighted_quantile(y(cell_j), w(cell_j), 0.5);
            else
                cell_val = sum(y(cell_j).*w(cell_j)) / wj;
            end
            num = num + pz(j) * cell_val;
            den = den + pz(j);
        end
    end
    if den > 0, yc(i) = num / den; end
end
good = ~isnan(xc) & ~isnan(yc);  xc = xc(good);  yc = yc(good);
end

function m = wmean(v, w)
% Weighted mean over strictly positive v only.
v = v(:); w = w(:);
ok = isfinite(v) & v > 0 & isfinite(w) & w > 0;
m = sum(v(ok).*w(ok)) / max(sum(w(ok)), 1e-300);
end

function draw_tilt(ax, xm, ys, keep, s_unc, XCUT, y_unc, y_end, CB, LW, FS)
% The used-capital value share against the shadow cost. Every coordinate is
% passed in; this function computes nothing, so the panel is identical whether
% it is drawn inside Figure 1 or on its own canvas.
CR = [0.55 0.55 0.55];
XL = -0.030 * XCUT;

plot(ax, [XL XCUT], [y_unc y_unc], ':',  'Color', CR, 'LineWidth', 0.9);
plot(ax, [0 0],     [0 100],       '--', 'Color', CR, 'LineWidth', 0.9);
plot(ax, [0 xm(keep)], 100*[s_unc ys(keep)], '-', 'Color', CB, 'LineWidth', LW);
plot(ax, 0, y_unc, 'o', 'MarkerEdgeColor', CB, 'MarkerFaceColor', 'w', ...
    'MarkerSize', 5.5, 'LineWidth', 1.3);

xa = XCUT*0.88;
plot(ax, [xa xa], [y_unc y_end], '-', 'Color', CB, 'LineWidth', 0.8);
plot(ax, xa, y_unc, 'v', 'MarkerEdgeColor', CB, 'MarkerFaceColor', CB, 'MarkerSize', 3.5);
plot(ax, xa, y_end, '^', 'MarkerEdgeColor', CB, 'MarkerFaceColor', CB, 'MarkerSize', 3.5);
text(ax, xa - 0.012*XCUT, 0.5*(y_unc + y_end), sprintf('+%.0f pp', y_end - y_unc), ...
    'Color', CB, 'FontSize', FS-1, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle');
text(ax, XCUT*0.035, y_unc - 0.6, 'unconstrained', 'Color', CR, ...
    'FontSize', FS-1, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax, XCUT*0.30, y_unc + 0.30*(y_end - y_unc), ...
    'tighter constraint \rightarrow more used capital', 'Color', CB, ...
    'FontSize', FS-1, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlabel(ax, 'Shadow cost of collateral constraint, \mu');
ylabel(ax, 'Used-capital value share (%)');
xlim(ax, [XL XCUT]);
ylim(ax, [y_unc - 0.18*(y_end - y_unc), y_end + 0.12*(y_end - y_unc)]);
grid(ax, 'on'); set(ax, 'GridAlpha', 0.10, 'Layer', 'bottom');
end
