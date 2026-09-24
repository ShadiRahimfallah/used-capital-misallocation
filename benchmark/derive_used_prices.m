function [qu, pu, Rused] = derive_used_prices(r, delta, kappa, Rnew, zeta, psi_u)
if nargin < 6 || isempty(psi_u), psi_u = 0; end

Rused = r + delta + kappa .* (1 - psi_u);

if isempty(zeta) || zeta <= 0
    qu = 0;
    pu = 0;
    return
end

qu = Rused .* (1 - (Rnew - r - delta) ./ zeta);

if qu <= 1e-12
    qu = 0;
    pu = 0;
else
    pu = qu ./ max(Rused, 1e-12);
end
end
