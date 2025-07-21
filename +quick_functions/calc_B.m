function B = calc_B(model, B_positions)
elem_centers = model.elem_centers;
elem_volumes = model.elem_volumes;

b = zeros(size(B_positions,1), 3, model.phantom.n_elec, 1);
progressbar
for ii = 1:model.phantom.n_elec
    elem_currents = model.elem_currents(:,:,ii);
    b(:,:,ii,1) = helpers.calc_B_at_points(B_positions, elem_centers, elem_currents, elem_volumes);
    progressbar(ii/model.phantom.n_elec)
end

b_cyl = zeros(size(b));
for ii = 1:model.phantom.n_elec
    b_cyl(:,:,ii,1) = quick_functions.from_cart_to_cyl(b(:,:,ii,1), B_positions);
end

bb = permute(b_cyl, [1,3,2]);
B = reshape(bb, 256, 3);


end