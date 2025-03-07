function [B] = calc_B_at_points(B_positions, elem_centers, e_curr, elem_vol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu0 = 4*pi*1e-7;
n_points = size(B_positions,1);
B = zeros(n_points, 3);

% precalculate current density * element volume
j_times_vol = e_curr.*elem_vol;

for i_point = 1:n_points
    point = B_positions(i_point,:);
    r = point - elem_centers;
    
    % if need to check for zero division
    r_mag = vecnorm(r.',2,1).';
    l_non_zero = find(r_mag);

    % dB = cross(j_times_vol(l_non_zero,:), r(l_non_zero,:))./norm(r(l_non_zero,:).^3); % magnetic flux density, in T
    dB = cross(j_times_vol(l_non_zero,:), r(l_non_zero,:))./vecnorm(r(l_non_zero,:),2,2).^3; % magnetic flux density, in T

    dB(isnan(dB)) = 0;
    B(i_point, :) = mu0/(4*pi)*sum(dB,1);
end

end