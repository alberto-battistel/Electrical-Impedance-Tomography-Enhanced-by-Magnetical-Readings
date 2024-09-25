home
clear
close all
init_eidors()


%%
n_elec = 16;
target.center = @(z) [0.0, 0.0, z];
target.radius = 0.02;
elec_vert_position = 0.1;
phantom_radius = 0.1;
phantom_height = 2*phantom_radius;
maxsz = 0.01;

el_pos = [-22.5, elec_vert_position; ...
        22.5, elec_vert_position; ...
        180-22.5, elec_vert_position; ...
        180+22.5, elec_vert_position];
el_sz  = [0.01,0,0.05].*ones(size(el_pos,1),3);

fmdl= ng_mk_cyl_models([phantom_height,phantom_radius,maxsz],el_pos,el_sz);
imdl = mk_common_model('a2c2',16); % Will replace most 
imdl.fwd_model = fmdl;
imdl.normalize_measurements = 0;
stim_pattern = mk_stim_patterns(size(el_pos,1),1,'{ad}','{ad}',{},10e-3);
imdl.fwd_model.stimulation = stim_pattern;


img_h = mk_image(imdl, 0.503); % muscle cond 1 MHz
img_h.fwd_solve.get_all_meas = 1;

show_fem(fmdl,[0,1.012])

%% positions where to calculate the magnetic field
n_points = 50;

[xx, yy, zz] = meshgrid(linspace(-phantom_radius, phantom_radius, n_points), phantom_radius, linspace(0, phantom_height, n_points));
B_positions = [xx(:), yy(:), zz(:)];

%
elem_centers = interp_mesh(imdl.fwd_model, 0); % center of elements
elem_areas = helpers.calc_area(imdl.fwd_model.nodes, imdl.fwd_model.elems);

vh = fwd_solve(img_h);
% vh = vh_.meas;
% bh_ = helpers.magnetic_acquisition(img_h, vh_, B_positions, elem_centers, elem_areas);
mu0 = 4*pi*1e-7;

n_measurements = length(vh.meas);
B = zeros(n_measurements, n_measurements, 3);
for i_volt = 1:n_measurements
    e_curr = calc_elem_current(img_h, vh.volt(:,i_volt));

    for i_point = 1:length(B_positions)
        point = B_positions(i_point,:);
        r = point - elem_centers;
        % if need to check for zero division
        r_mag = vecnorm(r.',2,1).';
        l = find(r_mag);
        
        dB = mu0/(4*pi)*cross(e_curr.*elem_areas, r)./norm(r.^3); % magnetic flux density, in T
        % dB = mu0/(4*pi)*cross_prod(e_curr(l).*elem_areas(l), r(l))./norm(r(l).^3); % magnetic flux density, in T
    
        B(i_point, i_volt, :) = sum(dB,1);
    end
end
% B = reshape(B,n_measurements*n_measurements, []);

%%
abs_B = vecnorm(B,2,3);

[xx, zz] = meshgrid(linspace(-phantom_radius, phantom_radius, n_points), linspace(0, phantom_height, n_points));

surf(xx,zz,reshape(abs_B(:,1),n_points,n_points))