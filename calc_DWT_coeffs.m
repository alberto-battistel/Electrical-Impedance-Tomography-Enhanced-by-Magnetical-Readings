function [coeffs_matrix] = calc_DWT_coeffs(n)
coeffs_matrix = zeros(prod(n),3);

idx = 0;
for ix = 1:2*n(1)
    for iy = 1:2*n(2)
        for iz = 1:2*n(3)
            idx = idx+1;
            coeffs_matrix(idx,:) = [ix, iy, iz];
        end
    end
end
end