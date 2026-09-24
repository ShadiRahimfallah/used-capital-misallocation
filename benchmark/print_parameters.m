function print_parameters(x, x_prev, scale)
names = {'W','r','Rnew','entry_cost','theta','se','rho','kappa','pw','sigma_y'};
fmts  = {'%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f'};

n = min(numel(x), numel(names));

have_prev = (nargin >= 3) && ~isempty(x_prev) && all(isfinite(x_prev));

cols_per_line = 3;
for i = 1:n
    if mod(i-1, cols_per_line) == 0
        fprintf('    ');
    end

    val_str = sprintf(fmts{i}, x(i));
    fprintf('%-6s = %s', names{i}, val_str);

    if have_prev && numel(x_prev) >= i && numel(scale) >= i
        dx      = x(i) - x_prev(i);
        dx_norm = abs(dx) / max(scale(i), 1e-12);
        if dx_norm > 1e-6
            fprintf(' (%+.2e)', dx);
        end
    end

    if mod(i, cols_per_line) == 0 || i == n
        fprintf('\n');
    else
        fprintf('   ');
    end
end
end
