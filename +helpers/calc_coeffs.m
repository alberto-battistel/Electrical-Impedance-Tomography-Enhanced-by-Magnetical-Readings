function [coeffs_matrix] = calc_coeffs(max_radial_coeff, max_cheb_coeff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

zern_coeffs = helpers.calc_zern_coeffs(max_radial_coeff);

n_zern_coeffs = size(zern_coeffs,1);

coeffs_matrix = zeros(n_zern_coeffs*(max_cheb_coeff+1),3);

coeffs_matrix(:,1:2) = repmat(zern_coeffs, max_cheb_coeff+1,1);

coeffs_matrix(:,3) = reshape(repmat(0:max_cheb_coeff,n_zern_coeffs,1),[],1);


end