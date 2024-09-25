function B = magnetic_acquisition(img, vh, coil_positions, elem_centers, elem_areas)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

mu0 = 4*pi*1e-7;
n_measurements = size(coil_positions,1);
if size(elem_centers,2) == 2
    cross_prod = @cross_prod_2d;
    B_dim = 1;
    B = zeros(n_measurements, n_measurements);
elseif size(elem_centers,2) == 3
    cross_prod = @cross;
    B_dim = 3;
    B = zeros(n_measurements, n_measurements, 3);
end

for i_volt = 1:n_measurements
    e_curr = calc_elem_current(img, vh.volt(:,i_volt));

    for i_point = 1:length(coil_positions)
        point = coil_positions(i_point,:);
        r = point - elem_centers;
        % if need to check for zero division
        % r_mag = norm(r);
    
        dB = mu0/(4*pi)*cross_prod(e_curr.*elem_areas, r)./norm(r.^3); % magnetic flux density, in T
    
        B(i_point, i_volt, :) = sum(dB,1);
    end
end
B = reshape(B,n_measurements*n_measurements, []);
end
