function [B] = calc_B_at_points(B_positions, elem_centers, e_curr, elem_vol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu0 = 4*pi*1e-7;
B = zeros(length(B_positions), 3);

for i_point = 1:length(B_positions)
    point = B_positions(i_point,:);
    r = point - elem_centers;
    
    % if need to check for zero division
    r_mag = vecnorm(r.',2,1).';
    l_non_zero = find(r_mag);

    dB = cross(e_curr(l_non_zero,:).*elem_vol(l_non_zero,:), r(l_non_zero,:))./norm(r(l_non_zero,:).^3); % magnetic flux density, in T

    B(i_point, :) = mu0/(4*pi)*sum(dB,1);
end

end