global theta eps_model egrid alpha eta W phi gamma Rnew qu ...
    productivity_weight F pu lambda_u

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
if ~isscalar(productivity_weight), productivity_weight = productivity_weight(1); end

lower_bound = max(eps_model^(1/4), 1e-3);
upper_bound = 10;

A = nodeunif(10000, lower_bound, upper_bound);
A = A(:);
N = numel(A);

Eeff = productivity_weight .* exp(egrid(1)) * ones(N,1);

[Kn, Ku, Kagg, L, Y, D] = compute_capital_allocation( ...
    A, Eeff, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, [], lambda_u);

valid   = isfinite(D);
A_valid = A(valid);
D_valid = D(valid);

if isempty(A_valid)
    Aumin = lower_bound;
else
    entry_condition = (D_valid >= W);

    if any(entry_condition)
        Aumin = min(A_valid(entry_condition));
    else
        Aumin = lower_bound;
    end
end
