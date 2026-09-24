global alpha eta W phi gamma Rnew qu delta r theta eps_model egrid ...
    productivity_weight kappa F pu lambda_u cf_cf1_zero_used_prices

A    = s(:,1);
Eeff = productivity_weight .* exp(egrid(s(:,2)));

[Kn, Ku, Kagg, L, Y, D, ~, binding] = compute_capital_allocation( ...
    A, Eeff, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, [], lambda_u);

if (~isempty(cf_cf1_zero_used_prices) && logical(cf_cf1_zero_used_prices)) ...
        || (abs(qu) <= 1e-14 && abs(pu) <= 1e-14)
    rho_ces       = 0;
    Psi_ces       = 1;
    leverage_gain = 1;
    R_bundle      = Rnew;
else
    rho_ces       = ((1 - phi) / phi) .* (Rnew ./ qu) .^ gamma;
    Psi_ces       = (phi^(1/gamma) + (1-phi)^(1/gamma) .* rho_ces.^((gamma-1)/gamma)) .^ (gamma/(gamma-1));
    leverage_gain = Psi_ces ./ (1 + rho_ces);
    R_bundle      = (Rnew + rho_ces .* qu) ./ (1 + rho_ces);
end
collateral = theta .* A;
