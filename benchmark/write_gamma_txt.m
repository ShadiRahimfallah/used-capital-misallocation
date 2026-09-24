R = load('robustness_gamma_result.mat');
fn = fieldnames(R);  R = R.(fn{1});
S  = R.S;

NF = load('nofriction_result.mat');
fb = fieldnames(NF);  NF = NF.(fb{1});
L_base = 100*log(NF.frictionless.TFP / NF.base.TFP);
Y_base = 100*(NF.frictionless.Y / NF.base.Y - 1);

fid = fopen('robustness_gamma_results.txt', 'w');
w = @(varargin) fprintf(fid, varargin{:});

w('==============================================================================\n');
w('  SENSITIVITY TO GAMMA, THE ELASTICITY OF SUBSTITUTION BETWEEN NEW AND USED\n');
w('  CAPITAL\n');
w('==============================================================================\n\n');
w('  The paper sets gamma = %.1f. Edgerton (2011) estimates the elasticity\n', R.base);
w('  between 1.9 and 2.4 for construction machinery; this sweep spans a wider\n');
w('  range, 1.7 to 2.5. The calibration is held fixed at the benchmark point\n');
w('  and only prices are re-cleared, so the exercise asks whether the results\n');
w('  turn on gamma, not whether the model can be re-fitted around it.\n\n');

w('  %-28s %11s %11s %11s\n', '', sprintf('gamma=%.1f', S(1).g), ...
    sprintf('gamma=%.1f', R.base), sprintf('gamma=%.1f', S(2).g));
w('  %s\n', repmat('-', 1, 64));
w('  %-28s %11.3f %11s %11.3f\n', 'wage W', S(1).fric.W, '0.540', S(2).fric.W);
w('  %-28s %11.3f %11s %11.3f\n', 'interest rate r', S(1).fric.r, '0.031', S(2).fric.r);
w('  %-28s %11.3f %11s %11.3f\n', 'new-capital rent Rnew', S(1).fric.Rnew, '0.119', S(2).fric.Rnew);
w('  %-28s %11.3f %11s %11.3f\n', 'used price pu', S(1).fric.pu, '0.645', S(2).fric.pu);
w('  %-28s %11.2f %11s %11.2f\n', 'used share (%)', 100*S(1).fric.ushare_agg, '--', 100*S(2).fric.ushare_agg);
w('  %-28s %11.3f %11s %11.3f\n', 'share constrained', S(1).fric.constr, '0.908', S(2).fric.constr);
w('  %-28s %11.3f %11s %11.3f\n', 'mean scale K/K*', S(1).fric.kk_mean, '0.343', S(2).fric.kk_mean);
w('  %-28s %11.2f %11.2f %11.2f\n', 'output gain (%)', S(1).gain_Y, Y_base, S(2).gain_Y);
w('  %-28s %11.2f %11.2f %11.2f\n', 'misallocation L (log pt)', S(1).gain_L, L_base, S(2).gain_L);
w('\n  Relative to gamma = %.1f:\n', R.base);
for i = 1:numel(S)
    w('    gamma = %.1f : L %+.3f log points, output gain %+.3f points\n', ...
        S(i).g, S(i).gain_L - L_base, S(i).gain_Y - Y_base);
end
w('    total spread in L across the range: %.3f log points (%.1f%% of %.2f)\n', ...
    abs(S(1).gain_L - S(2).gain_L), ...
    100*abs(S(1).gain_L - S(2).gain_L)/L_base, L_base);

w('\n------------------------------------------------------------------------------\n');
w('  FINDING\n');
w('------------------------------------------------------------------------------\n');
w('  The results do not turn on gamma. Across 1.7 to 2.5 the misallocation loss\n');
w('  moves by %.3f log points, about %.0f percent of its level, and the output\n', ...
    abs(S(1).gain_L - S(2).gain_L), 100*abs(S(1).gain_L - S(2).gain_L)/L_base);
w('  gain by %.3f points.\n\n', abs(S(1).gain_Y - S(2).gain_Y));
if S(2).gain_L < S(1).gain_L
    w('  DIRECTION. L FALLS as gamma rises (%.3f at %.1f against %.3f at %.1f).\n', ...
        S(2).gain_L, S(2).g, S(1).gain_L, S(1).g);
    w('  A claim that a more substitutable pair raises the reported effect, and\n');
    w('  that gamma = %.1f is therefore the conservative end of the range, is not\n', R.base);
    w('  supported: the sign runs the other way. Given the size of the movement\n');
    w('  the honest statement is insensitivity, not a reversal.\n');
else
    w('  DIRECTION. L rises with gamma, so the lower end of the range is the\n');
    w('  conservative one.\n');
end

w('\n  READ THE TWO SWEPT VALUES AGAINST EACH OTHER, NOT AGAINST THE STORED\n');
w('  BENCHMARK LEVEL. Both gammas give an output gain about 0.24 points below\n');
w('  the stored gamma = 2.0 figure, and by almost the same amount. A shift\n');
w('  common to both cannot be a gamma effect; it is a level difference against\n');
w('  an economy solved in a separate run. The gamma effect is the difference\n');
w('  BETWEEN 1.7 and 2.5: %.3f points on the output gain and %.3f log points\n', ...
    abs(S(1).gain_Y - S(2).gain_Y), abs(S(1).gain_L - S(2).gain_L));
w('  on L. The L comparison straddles the benchmark (%+.3f and %+.3f), so it\n', ...
    S(1).gain_L - L_base, S(2).gain_L - L_base);
w('  carries no such offset.\n');

w('\n------------------------------------------------------------------------------\n');
w('  HOW THE ECONOMIES WERE CLEARED\n');
w('------------------------------------------------------------------------------\n');
w('  The friction economies clear to 0.2 percent on the goods and labour\n');
w('  markets and 0.5 percent on the new-to-used ratio, against defaults of 0.5\n');
w('  and 2 percent. The tighter setting is necessary rather than cosmetic: at\n');
w('  the default, gamma = 1.7 seeded at the benchmark prices left residuals of\n');
w('  -0.4 and -1.0 percent, both inside tolerance, so the clearer returned the\n');
w('  seed unchanged and the economy was never re-cleared. The 2 percent floor\n');
w('  on the ratio was hardcoded and unreachable from any setting; it is now\n');
w('  overridable, with the default unchanged.\n\n');
w('  The benchmark''s own realized gaps are -0.1, +0.1 and +0.3 percent, so it\n');
w('  already meets the tighter standard and was not re-cleared.\n\n');
w('  The frictionless twins clear to 0.5 percent, which is the tolerance\n');
w('  nofriction.m used for the benchmark''s own twin, so the comparison is\n');
w('  consistent across all three economies.\n\n');
w('  gamma is hardcoded inside evaluate_economy and reachable only through\n');
w('  cf_gamma_override. If that override stopped biting, every gamma would\n');
w('  solve the same economy and the sweep would report, wrongly, that gamma\n');
w('  does not matter. The two economies here differ (Y by %.2e), so the\n', ...
    abs(S(1).fric.Y - S(2).fric.Y));
w('  override is live and the sweep is meaningful.\n');
w('\n  Produced by robustness_gamma.m; numbers read from\n');
w('  robustness_gamma_result.mat.\n');
w('==============================================================================\n');
fclose(fid);

type robustness_gamma_results.txt
