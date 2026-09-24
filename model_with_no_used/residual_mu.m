function f = residual_mu(mu, Eeff, coll, ...
    Rnew, qu, pu, phi, gamma, ...
    alpha, eta, W, lambda_u)
if nargin < 12 || isempty(lambda_u), lambda_u = 1.0; end

A_tfp = Eeff;

C_scalar = (alpha / W)^((alpha*eta)/(1-eta)) ...
    * (eta)^(1/(1-eta)) ...
    * (1-alpha)^((1-alpha*eta)/(1-eta));

Rn_tilde = Rnew + mu;
qu_tilde = qu   + mu * lambda_u * pu;

if abs(qu) <= 1e-14 && abs(pu) <= 1e-14
    error('residual_mu:noUsedCapital', ...
        ['residual_mu was called with qu = pu = 0. The single-vintage ' ...
         'economy is handled by compute_capital_allocation_kn_only.']);
end

rho_mu = ((1-phi)/phi) * (Rn_tilde / qu_tilde)^gamma;
rho_mu = max(rho_mu, 1e-12);

Psi = ( phi^(1/gamma) ...
    + (1-phi)^(1/gamma) * rho_mu^((gamma-1)/gamma) )^(gamma/(gamma-1));
if ~isfinite(Psi) || Psi <= 0, f = 1e10; return; end

R_eff = Rn_tilde + rho_mu * qu_tilde;
if R_eff <= 0, f = 1e10; return; end

pow1 = ((1-alpha)*eta) / (1-eta);
pow2 = -(1-alpha*eta)  / (1-eta);
Kn   = C_scalar * (Psi^pow1) * (R_eff^pow2) * Eeff^(1/(1-eta));
Kn   = max(Kn, 1e-12);

f = Kn * (1 + lambda_u * pu * rho_mu) - coll;
if ~isfinite(f), f = 1e10; end
end
