function [out1, out2] = menufun(flag, s, x)
global alpha eta W r egrid smine sminw smaxe smaxw theta F phi gamma Rnew qu ...
    eps_model pu ygrid Ey sigma_crra lambda_u productivity_weight

if isempty(sigma_crra), sigma_crra = 1.5; end
if isempty(productivity_weight), productivity_weight = 1; end

n = size(s,1);
A = s(:,1);

switch flag
    case 'be'
        Eeff = productivity_weight .* exp(egrid(s(:,2)));

        D = cca_profit_cached(A, Eeff, Rnew, qu, phi, gamma, ...
            alpha, eta, W, F, pu, theta, lambda_u);

        out1 = smine(1)*ones(n,1);
        out2 = min(D + (1+r)*A - eps_model^(1/3), smaxe(1));
    case 'bw'
        Dw   = worker_wage(s, W, ygrid, Ey, n);
        out1 = sminw(1)*ones(n,1);
        out2 = min(Dw + (1+r)*A - eps_model^(1/3), smaxw(1));
    case 'fe'
        Eeff = productivity_weight .* exp(egrid(s(:,2)));

        D = cca_profit_cached(A, Eeff, Rnew, qu, phi, gamma, ...
            alpha, eta, W, F, pu, theta, lambda_u);

        C    = D + (1+r)*A - x;
        C    = max(C, 1e-12);
        out1 = crra_utility(C, sigma_crra);
        out2 = [];
    case 'fw'
        Dw   = worker_wage(s, W, ygrid, Ey, n);
        C    = Dw + (1+r)*A - x;
        C    = max(C, 1e-12);
        out1 = crra_utility(C, sigma_crra);
        out2 = [];
    otherwise
        error('menufun: unknown flag');
end
end

function D = cca_profit_cached(A, Eeff, Rnew, qu, phi, gamma, ...
    alpha, eta, W, F, pu, theta, lambda_u)
global MENU_CCA_CACHE_OFF cf_cf1_zero_used_prices

persistent key_A key_E key_sig D_cache

cf1 = (~isempty(cf_cf1_zero_used_prices) && logical(cf_cf1_zero_used_prices));
sig = [Rnew; qu; phi; gamma; alpha; eta; W; F; pu; theta; lambda_u; double(cf1)];

use_cache = isempty(MENU_CCA_CACHE_OFF) || ~MENU_CCA_CACHE_OFF;

if use_cache && ~isempty(D_cache) ...
        && numel(D_cache) == numel(A) ...
        && isequal(key_sig, sig) ...
        && isequal(key_A, A) ...
        && isequal(key_E, Eeff)
    D = D_cache;
    return;
end

[~, ~, ~, ~, ~, D] = compute_capital_allocation( ...
    A, Eeff, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, ...
    'vectorized', lambda_u);

if use_cache
    key_A   = A;
    key_E   = Eeff;
    key_sig = sig;
    D_cache = D;
end
end

function Dw = worker_wage(s, W, ygrid, Ey, n)
if size(s,2) >= 3 && ~isempty(ygrid)
    iy = s(:,3);
    Dw = W .* ygrid(iy);
else
    Dw = W * ones(n, 1);
end
end

function u = crra_utility(c, sigma)
if abs(sigma - 1) < 1e-8
    u = log(c);
else
    u = (c.^(1 - sigma) - 1) / (1 - sigma);
end
end
