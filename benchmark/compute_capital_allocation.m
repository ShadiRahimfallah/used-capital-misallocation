function [Kn, Ku, Kagg, L, Y, D, mu, binding] = compute_capital_allocation( ...
    A, Eeff, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, method, lambda_u)
global cf_cf1_zero_used_prices

if nargin < 13 || isempty(method), method = 'serial'; end
if nargin < 14 || isempty(lambda_u), lambda_u = 1.0; end

N    = numel(A);
A    = A(:);
Eeff = Eeff(:);

if ~isscalar(theta),    theta    = theta(1);    end
if ~isscalar(lambda_u), lambda_u = lambda_u(1); end
if ~isscalar(pu),       pu       = pu(1);       end
if ~isscalar(phi),      phi      = phi(1);      end
if ~isscalar(gamma),    gamma    = gamma(1);    end
if ~isscalar(alpha),    alpha    = alpha(1);    end
if ~isscalar(eta),      eta      = eta(1);      end
if ~isscalar(W),        W        = W(1);        end
if ~isscalar(F),        F        = F(1);        end
if ~isscalar(Rnew),     Rnew     = Rnew(1);     end
if ~isscalar(qu),       qu       = qu(1);       end

cf1_kn_only = (~isempty(cf_cf1_zero_used_prices) && logical(cf_cf1_zero_used_prices)) ...
    || (abs(qu) <= 1e-14 && abs(pu) <= 1e-14);
if cf1_kn_only
    error('compute_capital_allocation:noUsedCapital', ...
        ['This is the benchmark economy, which always has a second-hand ' ...
         'market. For the single-vintage economy use model_with_no_used/.']);
end

qu_alg  = qu;
rho_ces = ((1 - phi) / phi) .* (Rnew ./ qu_alg) .^ gamma;
Psi_ces = ( phi^(1/gamma) ...
    + (1-phi)^(1/gamma) .* rho_ces.^((gamma-1)/gamma) ) ...
    .^ (gamma/(gamma-1));
R_eff = Rnew + rho_ces .* qu_alg;

C_prod = (alpha ./ W) .^ ((alpha*eta)/(1-eta)) ...
    .* (eta) .^ (1/(1-eta)) ...
    .* (1-alpha) .^ ((1-alpha*eta)/(1-eta));

pow1 = ((1-alpha)*eta) / (1-eta);
pow2 = -(1-alpha*eta) / (1-eta);

Kn_uc   = C_prod .* (Psi_ces .^ pow1) .* (R_eff .^ pow2) .* Eeff.^(1/(1-eta));
Ku_uc   = rho_ces .* Kn_uc;
Kagg_uc = Psi_ces .* Kn_uc;

L_uc = ( (alpha*eta./W) .* Eeff ...
    .* (Psi_ces .* Kn_uc).^((1-alpha)*eta) ).^(1/(1-alpha*eta));

Y_uc = Eeff .* (L_uc.^alpha .* Kagg_uc.^(1-alpha)).^eta;
D_uc = Y_uc - W.*L_uc - Rnew.*Kn_uc - qu.*Ku_uc - W.*F;

collateral = theta .* A;
binding    = (Kn_uc + lambda_u .* pu .* Ku_uc) > collateral;

Kn   = Kn_uc;
Ku   = Ku_uc;
Kagg = Kagg_uc;
L    = L_uc;
Y    = Y_uc;
D    = D_uc;
mu   = zeros(N, 1);

if ~any(binding), return; end

bind_idx = find(binding);
n_bind   = numel(bind_idx);

switch lower(method)
    case 'vectorized'
        Eeff_b = Eeff(bind_idx);
        coll_b = collateral(bind_idx);

        mu_lo = zeros(n_bind, 1);
        mu_hi = 500 * ones(n_bind, 1);

        for iter_bi = 1:60
            mu_mid = (mu_lo + mu_hi) / 2;

            den_c = max(qu + mu_mid .* lambda_u .* pu, 1e-14);
            rho_c = ((1-phi)/phi) .* ((Rnew + mu_mid) ./ den_c).^gamma;
            Psi_c = ( phi^(1/gamma) ...
                + (1-phi)^(1/gamma) .* rho_c.^((gamma-1)/gamma) ) ...
                .^(gamma/(gamma-1));
            R_eff_c = (Rnew + mu_mid) + rho_c .* den_c;

            Kn_foc = C_prod .* Psi_c.^pow1 .* R_eff_c.^pow2 .* Eeff_b.^(1/(1-eta));

            g = Kn_foc .* (1 + lambda_u .* pu .* rho_c) - coll_b;

            too_high         = (g > 0);
            mu_lo(too_high)  = mu_mid(too_high);
            mu_hi(~too_high) = mu_mid(~too_high);
        end

        mu_sol = mu_mid;

        den_cf = max(qu + mu_sol .* lambda_u .* pu, 1e-14);
        rho_cf = ((1-phi)/phi) .* ((Rnew + mu_sol) ./ den_cf).^gamma;
        Psi_cf = ( phi^(1/gamma) ...
            + (1-phi)^(1/gamma) .* rho_cf.^((gamma-1)/gamma) ) ...
            .^(gamma/(gamma-1));

        Kn_c   = coll_b ./ (1 + lambda_u .* pu .* rho_cf);
        Ku_c   = rho_cf .* Kn_c;
        Kagg_c = Psi_cf .* Kn_c;

        A_tfp_b = Eeff_b;
        L_c = ( (alpha*eta/W) .* A_tfp_b ...
            .* (Psi_cf .* Kn_c).^((1-alpha)*eta) ).^(1/(1-alpha*eta));
        Y_c = A_tfp_b .* (L_c.^alpha .* Kagg_c.^(1-alpha)).^eta;
        D_c = Y_c - W.*L_c - Rnew.*Kn_c - qu.*Ku_c - W.*F;

        ok_b = isfinite(mu_sol) & isfinite(Kn_c) & (Kn_c > 0) ...
            & isfinite(Y_c) & isfinite(D_c);

        g_final = Kn_c .* (1 + lambda_u .* pu .* rho_cf) - coll_b;
        ok_b    = ok_b & (abs(g_final) <= 1e-6 * max(coll_b, 1));
        if any(~ok_b)
            warning('compute_capital_allocation:muCeiling', ...
                ['%d of %d constrained firms did not satisfy the collateral ' ...
                 'constraint after bisection. Raise the mu ceiling above 500 ' ...
                 'in this function and re-solve.'], sum(~ok_b), numel(ok_b));
        end

        Kn(bind_idx(ok_b))   = Kn_c(ok_b);
        Ku(bind_idx(ok_b))   = Ku_c(ok_b);
        Kagg(bind_idx(ok_b)) = Kagg_c(ok_b);
        L(bind_idx(ok_b))    = L_c(ok_b);
        Y(bind_idx(ok_b))    = Y_c(ok_b);
        D(bind_idx(ok_b))    = D_c(ok_b);
        mu(bind_idx(ok_b))   = mu_sol(ok_b);
    case 'parfor'
        Kn_b   = zeros(n_bind,1);
        Ku_b   = zeros(n_bind,1);
        Kagg_b = zeros(n_bind,1);
        L_b    = zeros(n_bind,1);
        Y_b    = zeros(n_bind,1);
        D_b    = zeros(n_bind,1);
        mu_b   = nan(n_bind,1);

        parfor ii = 1:n_bind
            [Kn_b(ii), Ku_b(ii), Kagg_b(ii), L_b(ii), ...
                Y_b(ii), D_b(ii), mu_b(ii)] = solve_one_constrained( ...
                Eeff(bind_idx(ii)), collateral(bind_idx(ii)), ...
                Rnew, qu, pu, phi, gamma, alpha, eta, W, F, lambda_u);
        end

        ok = isfinite(mu_b);
        Kn(bind_idx(ok))   = Kn_b(ok);
        Ku(bind_idx(ok))   = Ku_b(ok);
        Kagg(bind_idx(ok)) = Kagg_b(ok);
        L(bind_idx(ok))    = L_b(ok);
        Y(bind_idx(ok))    = Y_b(ok);
        D(bind_idx(ok))    = D_b(ok);
        mu(bind_idx(ok))   = mu_b(ok);
    otherwise
        for ii = 1:n_bind
            idx = bind_idx(ii);
            [Kn_c, Ku_c, Kagg_c, L_c, Y_c, D_c, mu_sol] = ...
                solve_one_constrained( ...
                Eeff(idx), collateral(idx), ...
                Rnew, qu, pu, phi, gamma, alpha, eta, W, F, lambda_u);
            if ~isfinite(mu_sol), continue; end
            Kn(idx)   = Kn_c;
            Ku(idx)   = Ku_c;
            Kagg(idx) = Kagg_c;
            L(idx)    = L_c;
            Y(idx)    = Y_c;
            D(idx)    = D_c;
            mu(idx)   = mu_sol;
        end
end
end

function [Kn_c, Ku_c, Kagg_c, L_c, Y_c, D_c, mu_out] = ...
    solve_one_constrained(Eeff_i, coll_i, ...
    Rnew, qu, pu, phi, gamma, alpha, eta, W, F, lambda_u)
if nargin < 12 || isempty(lambda_u), lambda_u = 1.0; end
Kn_c   = 0;
Ku_c   = 0;
Kagg_c = 0;
L_c    = 0;
Y_c    = 0;
D_c    = 0;
mu_out = NaN;

mu_sol = solve_mu_bisect(Eeff_i, coll_i, ...
    Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);
if ~isfinite(mu_sol) || mu_sol < 0, return; end

den_1 = max(qu + mu_sol*lambda_u*pu, 1e-14);
rho_c = ((1-phi)/phi) * ((Rnew + mu_sol) / den_1)^gamma;
if ~isfinite(rho_c) || rho_c < 0, return; end

Psi_c = ( phi^(1/gamma) ...
    + (1-phi)^(1/gamma) * rho_c^((gamma-1)/gamma) )^(gamma/(gamma-1));
if ~isfinite(Psi_c) || Psi_c <= 0, return; end

Kn_c   = coll_i / (1 + lambda_u * pu * rho_c);
Ku_c   = rho_c * Kn_c;
Kagg_c = Psi_c * Kn_c;
if ~(isfinite(Kn_c) && Kn_c > 0 && isfinite(Kagg_c) && Kagg_c > 0)
    mu_out = NaN;
    return;
end

A_tfp_i = Eeff_i;
L_c = ((alpha*eta / W) * A_tfp_i ...
    * (Psi_c * Kn_c)^((1-alpha)*eta))^(1/(1-alpha*eta));
if ~isfinite(L_c) || L_c <= 0, mu_out = NaN; return; end

Y_c = A_tfp_i * (L_c^alpha * Kagg_c^(1-alpha))^eta;
D_c = Y_c - W*L_c - Rnew*Kn_c - qu*Ku_c - W*F;
if ~(isfinite(Y_c) && isfinite(D_c)), mu_out = NaN; return; end

mu_out = mu_sol;
end
