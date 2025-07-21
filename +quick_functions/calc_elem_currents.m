function elem_currents = calc_elem_currents(img, vh)
n_stimulations = 16;
elem_currents = zeros(length(img.elem_data),3,n_stimulations);
for ii = 1:n_stimulations
    elem_currents(:,:,ii) = calc_elem_current(img, vh.volt(:,ii));
end
end