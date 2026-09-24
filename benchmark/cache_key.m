function key = cache_key(opts)
if nargin < 1, opts = struct(); end
gv = @(f,d) getdef_local(opts, f, d);

nm = {'theta','eta_p','omega','kappa','sigma_y','F','rho_y','beta','lambda_u'};
dv = [1.97, 5.70, 0.855, 0.100, 0.275, 0.0, 0.95, 0.865, 1.0];

parts = cell(1, numel(nm));
for j = 1:numel(nm)
    parts{j} = sprintf('%s=%.10g', nm{j}, gv(nm{j}, dv(j)));
end

key = sprintf('%s|fast=%d|clear=%d', strjoin(parts, ' '), ...
    logical(gv('fast', true)), logical(gv('do_clear', true)));
end

function v = getdef_local(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
