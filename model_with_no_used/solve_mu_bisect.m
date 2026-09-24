function mu_sol = solve_mu_bisect(Eeff, coll, ...
    Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u)
if nargin < 11 || isempty(lambda_u), lambda_u = 1.0; end

tol    = 1e-11;
maxexp = 60;
maxbis = 200;

a_lo = 0;
b_hi = 1.0;
fa   = residual_mu(a_lo, Eeff, coll, Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);
if ~isfinite(fa), fa = 1; end

fb = residual_mu(b_hi, Eeff, coll, Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);

it = 0;
while isfinite(fb) && sign(fb) == sign(fa) && b_hi < 1e6 && it < maxexp
    b_hi = 2 * b_hi;
    fb   = residual_mu(b_hi, Eeff, coll, Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);
    it   = it + 1;
end

if ~(isfinite(fb) && sign(fb) ~= sign(fa))
    mu_sol = NaN;
    return;
end

for iter = 1:maxbis
    if (b_hi - a_lo) <= tol
        break;
    end
    mid = 0.5 * (a_lo + b_hi);
    fx  = residual_mu(mid, Eeff, coll, Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);
    if fx > 0
        a_lo = mid;
    else
        b_hi = mid;
    end
end

mu_sol = max(0.5 * (a_lo + b_hi), 0);

f_final = residual_mu(mu_sol, Eeff, coll, Rnew, qu, pu, phi, gamma, alpha, eta, W, lambda_u);
if ~isfinite(f_final) || abs(f_final) > 1e-4 * max(coll, 1)
    mu_sol = NaN;
end
end
