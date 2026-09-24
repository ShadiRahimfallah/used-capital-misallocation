function [Kn, Ku, Kagg, L, Y, D, mu, binding] = compute_capital_allocation_kn_only( ...
    A, Eeff, Rnew, W, F, alpha, eta, theta, method, lambda_u, qu, pu)
N    = numel(A);
A    = A(:);
Eeff = Eeff(:);

C_prod = (alpha ./ W) .^ ((alpha*eta)/(1-eta)) ...
    .* (eta) .^ (1/(1-eta)) ...
    .* (1-alpha) .^ ((1-alpha*eta)/(1-eta));
pow2  = -(1-alpha*eta) / (1-eta);
Kn_uc = C_prod .* (Rnew .^ pow2) .* Eeff.^(1/(1-eta));

collateral = theta .* A;
binding    = Kn_uc > collateral;

Kn          = Kn_uc;
Kn(binding) = collateral(binding);
Kn          = max(Kn, 0);

L = ( (alpha*eta./W) .* Eeff .* Kn.^((1-alpha)*eta) ).^(1/(1-alpha*eta));
Y = Eeff .* (L.^alpha .* Kn.^(1-alpha)).^eta;
D = Y - W.*L - Rnew.*Kn - W.*F;

mu          = zeros(N, 1);
mu(binding) = max(0, (1-alpha).*eta.*Y(binding) ./ max(Kn(binding), 1e-14) - Rnew);

Ku   = zeros(N, 1);
Kagg = Kn;
end
