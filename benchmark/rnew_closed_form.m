function [Rnew, s, qu_, pu_] = rnew_closed_form(r)
global delta zeta kappa phi gamma KnKu_ratio psi_u

psi_loc = 0;
if ~isempty(psi_u), psi_loc = psi_u; end

KnKu_loc = 2.5699;
if ~isempty(KnKu_ratio), KnKu_loc = KnKu_ratio; end

if isempty(zeta) || zeta <= 0 || phi >= 1 - 1e-12
    s    = 0;
    Rnew = r + delta;
    qu_  = 0;
    pu_  = 0;
    return
end

Rused = r + delta + kappa * (1 - psi_loc);
m     = ( (phi/(1-phi)) / KnKu_loc ) ^ (1/gamma);

s = ( m*Rused - (r + delta) ) / ( 1 + m*Rused/zeta );
s = min(max(s, 1e-5), zeta - 1e-5);

Rnew = r + delta + s;
pu_  = 1 - s/zeta;
qu_  = Rused * pu_;
end
