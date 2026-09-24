function [x, xcoord] = nodeunif(n, a, b)
d = length(n);
if d == 1
    x      = linspace(a, b, n)';
    xcoord = x;
else
    xcoord = cell(1, d);
    for i = 1:d
        xcoord{i} = linspace(a(i), b(i), n(i))';
    end
    x = gridmake(xcoord{:});
end
end
