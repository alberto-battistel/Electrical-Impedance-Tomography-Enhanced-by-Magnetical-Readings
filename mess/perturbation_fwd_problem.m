function [J_e, J_m, A] = perturbation_fwd_problem(mdl, pert_amplitude, zern_order, coil_positions)
% given a model mdl calculate the discrete jacobian J and associated 
% reconstruction matrix A for a pert_amplitude
% with n_freqs and a tensor product of freqs_order and zern_order 

un_perturbed_value = 1;

elem_centers = interp_mesh(mdl.fwd_model, 0); % center of elements
cyl_elem_centers = cylindrical_elem_centers(elem_centers);

[coefficients_matrix] = calc_zern_coeffs(zern_order);
n_space_pert = size(coefficients_matrix,1);
space_basis_set = zernfun(coefficients_matrix(:,1), coefficients_matrix(:,2), cyl_elem_centers(:,1), cyl_elem_centers(:,2), 'norm');

A = space_basis_set;
n_combs = size(A,2);

%%

homo_img = mk_image(mdl, un_perturbed_value);
homo_img.fwd_solve.get_all_meas = 1;
% homo_img.show_slices.do_colourbar = true;
homo_elem_data = homo_img.elem_data;

homo_vh = fwd_solve(homo_img);
v0 = homo_vh.meas;

pert_img = cell(n_space_pert ,1);
fwd_pert = cell(n_space_pert ,1);

for i_combs = 1:n_combs
    aa_ = A(:,i_combs);

    pert_img{i_combs} = homo_img;
    pert_img{i_combs}.elem_data = homo_elem_data + pert_amplitude*aa_;
    fwd_pert{i_combs} = fwd_solve(pert_img{i_combs});
    
end
v_pert = reshape(cell2mat(cellfun(@(c) c.meas, fwd_pert, 'UniformOutput', false)), 208, []);


%%

v_diff = v_pert-v0;

J_e = v_diff/pert_amplitude;

%% magnetic part
elem_areas = calc_area(mdl.fwd_model.nodes, mdl.fwd_model.elems);

B0 = magnetic_acquisition(homo_img, homo_vh, coil_positions, elem_centers, elem_areas);
B = zeros(length(coil_positions).^2, n_combs);
for i_combs = 1:n_combs
    B(:,i_combs) = magnetic_acquisition(pert_img{i_combs}, fwd_pert{i_combs}, coil_positions, elem_centers, elem_areas);
end

B_diff = B-B0;
J_m = B_diff/pert_amplitude;

end


