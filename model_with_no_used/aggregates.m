function AG = aggregates()
global INV_ME INV_MW INV_AG egrid ygrid k Ny productivity_weight W r Rnew qu ...
    pu phi gamma alpha eta theta F lambda_u

if isempty(productivity_weight), productivity_weight = 1; end

ag  = INV_AG(:);
S   = gridmake(ag, (1:k)', (1:Ny)');
a_s = S(:,1);
z_s = S(:,2);
y_s = S(:,3);

ne  = INV_ME(:);
nw  = INV_MW(:);
tot = sum(ne) + sum(nw);
if tot > 0, ne = ne/tot;  nw = nw/tot; end
mass_e = sum(ne);

Eeff = productivity_weight .* exp(egrid(z_s));
[Kn_e, Ku_e, Kagg_e, L_e, Y_e, ~, mu_e, bind_e] = ...
    compute_capital_allocation(a_s, Eeff, Rnew, qu, phi, gamma, ...
    alpha, eta, W, F, pu, theta, 'vectorized', lambda_u);
Debt_e = max(0, Kn_e + pu .* Ku_e - a_s);

AG.Y    = ne' * Y_e;
AG.Kn   = ne' * Kn_e;
AG.Ku   = ne' * Ku_e;
AG.Kagg = ne' * Kagg_e;
AG.L    = ne' * L_e;
AG.Debt = ne' * Debt_e;
AG.Lsup = nw' * ygrid(y_s);
AG.ent  = mass_e;
AG.KnKu = AG.Kn / max(AG.Ku, 1e-12);
AG.DtoY = AG.Debt / max(AG.Y, 1e-12);

AG.TFP  = AG.Y / max((AG.L^alpha * AG.Kagg^(1-alpha))^eta, 1e-12);

w_e = ne / max(mass_e, 1e-12);
AG.constr  = sum(w_e(bind_e > 0));
AG.mu_mean = w_e' * mu_e;

scale_i = Eeff .^ (1/(1-eta));
unc     = (~bind_e) & (ne > 0);
if sum(ne(unc)) > 1e-8
    c_i   = Kn_e(unc) ./ max(scale_i(unc), 1e-12);
    qs    = quantile(c_i, [0.25 0.50 0.75]);
    c_hat = qs(2);
    AG.kk_c     = c_hat;
    AG.kk_c_iqr = qs(3) - qs(1);
    kk_i   = min(Kn_e ./ max(c_hat .* scale_i, 1e-12), 1.2);
    AG.kk_mean = w_e' * kk_i;
else
    AG.kk_c = NaN;
    AG.kk_c_iqr = NaN;
    AG.kk_mean = NaN;
end

AG.W = W;
AG.r = r;
AG.Rnew = Rnew;
AG.qu = qu;
AG.pu = pu;
AG.theta = theta;
end
