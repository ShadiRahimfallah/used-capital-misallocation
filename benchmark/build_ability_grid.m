function [egrid, pi_z, P, Perg, zP, k] = build_ability_grid(eta_p, omega, fast_diagnostic, legacy)
if nargin < 3 || isempty(fast_diagnostic), fast_diagnostic = false; end
if nargin < 4 || isempty(legacy),          legacy          = false; end

ap_safe = max(eta_p, 1.10);
om_safe = max(min(omega, 0.999), 0.0);

if legacy
    if fast_diagnostic, kk = 21; else, kk = 40; end
    F_all = [linspace(0.633, 0.998, kk-2)'; 0.999; 0.9995];
else
    if fast_diagnostic
        n_bot = 15;  n_tail = 15;
    else
        n_bot = 24;  n_tail = 24;
    end
    F_bot  = linspace(0.02,  0.61,  n_bot)';
    F_body = linspace(0.633, 0.998, n_tail)';
    F_all  = [F_bot; F_body; 0.999; 0.9995];
end

k = numel(F_all);

zP_raw  = (1 - F_all).^(-1/ap_safe);

F_edges = [0; 0.5*(F_all(1:end-1) + F_all(2:end)); 1];
pi_z    = diff(F_edges);
pi_z    = pi_z / sum(pi_z);

zP    = zP_raw / (pi_z' * zP_raw);
egrid = log(zP);

P    = (1 - om_safe) * repmat(pi_z', k, 1) + om_safe * eye(k);
Perg = pi_z;

assert(all(abs(sum(P,2) - 1) < 1e-10), 'build_ability_grid: P rows must sum to 1');
assert(abs(sum(Perg) - 1) < 1e-10,     'build_ability_grid: Perg must sum to 1');
assert(abs(pi_z' * zP - 1) < 1e-8,     'build_ability_grid: E[z]=1 normalisation failed');
end
