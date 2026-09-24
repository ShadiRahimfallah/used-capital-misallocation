function Mout = inv_shock(Min, Na, k, Ny, P, Py)
M3 = reshape(Min, [Na, k, Ny]);
for iy = 1:Ny
    M3(:,:,iy) = M3(:,:,iy) * P;
end
M2   = reshape(M3, [Na*k, Ny]) * Py;
Mout = reshape(M2, [Na*k*Ny, 1]);
end
