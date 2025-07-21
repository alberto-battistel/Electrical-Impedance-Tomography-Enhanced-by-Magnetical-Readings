function [obj] = make_model(phantom, current_ampl)

el_pos = [-360/phantom.n_elec/2+(0:phantom.n_elec-1).'/phantom.n_elec*360,phantom.elec_vert_position.*ones(phantom.n_elec,1)];
el_sz  = [phantom.elec_radius, 0, phantom.max_el_sz].*ones(size(el_pos,1),3);


if isfield(phantom, 'extra')
    fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz, phantom.extra);
else
    fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz);
end

phantom.with_extras = false;
if length(fwd_model.mat_idx) > 1
    phantom.with_extras = true;
end

obj.phantom = phantom;
obj.current_ampl = current_ampl;

stim_pattern = mk_stim_patterns(size(el_pos,1), 1, '{ad}', '{ad}', {}, current_ampl);
            
fwd_model.stimulation = stim_pattern;

img = mk_image(fwd_model, phantom.background);
img.show_slices.levels = [inf, inf, phantom.elec_vert_position];

if phantom.with_extras
    img.elem_data(fwd_model.mat_idx{2}) = phantom.ball;
end

img.fwd_solve.get_all_meas = 1;
obj.img = img;

obj.elem_centers = interp_mesh(fwd_model, 0); % center of elements
obj.elem_volumes = get_elem_volume(fwd_model);

end