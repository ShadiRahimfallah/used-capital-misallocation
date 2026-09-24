function [ygrid, Py, Ey, piy] = rouwenhorst_y(Ny, rho_y, sigma_y)
if nargin < 1 || isempty(Ny),      Ny = 7;       end
if nargin < 2 || isempty(rho_y),   rho_y = 0.95; end
if nargin < 3 || isempty(sigma_y), sigma_y = 0.30; end

p = (1 + rho_y) / 2;
P = [p, 1-p; 1-p, p];
for n = 3:Ny
    Pn               = zeros(n);
    Pn(1:n-1, 1:n-1) = Pn(1:n-1, 1:n-1) +     p   * P;
    Pn(1:n-1, 2:n  ) = Pn(1:n-1, 2:n  ) + (1-p)   * P;
    Pn(2:n,   1:n-1) = Pn(2:n,   1:n-1) + (1-p)   * P;
    Pn(2:n,   2:n  ) = Pn(2:n,   2:n  ) +     p   * P;
    Pn(2:n-1, :    ) = Pn(2:n-1, :    ) / 2;
    P                = Pn;
end
Py = P;

sigma_uncond = sigma_y / sqrt(1 - rho_y^2);
zmax         = sigma_uncond * sqrt(Ny - 1);
ygrid_log    = linspace(-zmax, zmax, Ny)';
ygrid        = exp(ygrid_log);

[V, D] = eig(Py');
[~, k] = min(abs(diag(D) - 1));
piy = abs(V(:, k));
piy = piy / sum(piy);

Ey = piy(:)' * ygrid(:);
end
