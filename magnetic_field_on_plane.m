home
clear
close all
init_eidors()


%%

phantom.n_elec = 16;
phantom.elec_vert_position = 0.1;
phantom.phantom_radius = 0.1;
phantom.phantom_height = 2*phantom.phantom_radius;

maxsz = 0.02;
max_el_sz = maxsz;

fprintf('maxsz %f\n', maxsz)
fprintf('max_el_sz %f\n', max_el_sz)

[fmdl, img_h, vh, elem_centers, elem_volumes, e_curr] = make_model_and_get_all(phantom, maxsz, max_el_sz);

strct = struct('B_positions_struct', [], 'B', [], 'abs_B', []);
B_results.vert_plane_tangential = strct;
B_results.horz_plane_tangential = strct;
B_results.vert_plane_sagital = strct;

n_points = 25;

what_2_plot.fem = true;
what_2_plot.abs_B = true;
what_2_plot.B = true;

%% vert plane tangential to phantom
B_positions_struct = calc_position_vertical_plane_tangential_to_phantom(n_points, phantom);
B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
abs_B = vecnorm(B,2,2);

make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot)

B_results.vert_plane_tangential.B_positions_struct = B_positions_struct;
B_results.vert_plane_tangential.B = B;
B_results.vert_plane_tangential.abs_B = abs_B;

%% horz plane tangential to phantom
B_positions_struct = calc_position_horz_plane_tangential_to_phantom(n_points, phantom);
B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
abs_B = vecnorm(B,2,2);

make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot)

B_results.horz_plane_tangential.B_positions_struct = B_positions_struct;
B_results.horz_plane_tangential.B = B;
B_results.horz_plane_tangential.abs_B = abs_B;

%% vert plane sagital to phantom
[B_positions_struct] = calc_position_vert_plane_sagital_to_phantom(n_points, phantom);
B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
abs_B = vecnorm(B,2,2);

make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot)

B_results.vert_plane_sagital.B_positions_struct = B_positions_struct;
B_results.vert_plane_sagital.B = B;
B_results.vert_plane_sagital.abs_B = abs_B;









%% function declarations

function [fmdl, img_h, vh, elem_centers, elem_volumes, e_curr] = make_model_and_get_all(phantom,maxsz,max_el_sz)
n_elec = phantom.n_elec;
elec_vert_position = phantom.elec_vert_position;
phantom_radius = phantom.phantom_radius;
phantom_height = phantom.phantom_height;

el_pos = [-360/n_elec/2+(0:n_elec-1).'/n_elec*360,elec_vert_position.*ones(16,1)];
el_sz  = [0.01,0,max_el_sz].*ones(size(el_pos,1),3);

fmdl = ng_mk_cyl_models([phantom_height,phantom_radius,maxsz], el_pos, el_sz);
imdl = mk_common_model('a2c2',16); % Will replace most 
imdl.fwd_model = fmdl;
imdl.normalize_measurements = 0;
stim_pattern = mk_stim_patterns(size(el_pos,1),1,'{ad}','{ad}',{},10e-3);

imdl.fwd_model.stimulation = stim_pattern;

img_h = mk_image(imdl, 0.503); % muscle cond 1 MHz
img_h.fwd_solve.get_all_meas = 1;

% figure()
% show_fem(fmdl,[0,1.012])

%%
%
elem_centers = interp_mesh(imdl.fwd_model, 0); % center of elements
elem_volumes = helpers.calc_element_volume(imdl.fwd_model.elems, imdl.fwd_model.nodes);

vh = fwd_solve(img_h);

% figure()
% show_current(img_h,vh.volt(:,1));

e_curr = calc_elem_current(img_h, vh.volt(:,1));

end

function B_position_plot3(B_positions)
plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), 'ro','MarkerFaceColor','r');
end

function [B_positions_struct] = calc_position_vertical_plane_tangential_to_phantom(n_points, phantom)
x_dir = linspace(-phantom.phantom_radius, phantom.phantom_radius, n_points);
y_dir = phantom.phantom_radius;
z_dir = linspace(0, phantom.phantom_height, n_points);
[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];

B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);

[xx, zz] = meshgrid(x_dir, z_dir);

B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
end

function [B_positions_struct] = calc_position_horz_plane_tangential_to_phantom(n_points, phantom)
x_dir = linspace(-phantom.phantom_radius, phantom.phantom_radius, n_points);
y_dir = linspace(phantom.phantom_radius, 3*phantom.phantom_radius, n_points);
z_dir = phantom.phantom_height/2;

[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];

B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);

[xx, zz] = meshgrid(x_dir, y_dir);

B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
end

function [B_positions_struct] = calc_position_vert_plane_sagital_to_phantom(n_points, phantom)
x_dir = 0;
y_dir = linspace(phantom.phantom_radius, 3*phantom.phantom_radius, n_points);
z_dir = linspace(0, phantom.phantom_height, n_points);

[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];

B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);

[xx, zz] = meshgrid(z_dir, y_dir);

B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
end

function make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot)

if what_2_plot.fem
    figure()
    clf
    hold on
    show_fem(fmdl,[0,1.012])
    B_positions_struct.plot3()
    hold off
end

if what_2_plot.abs_B
    figure()
    B_positions_struct.surf(abs_B)
end

if what_2_plot.B
    figure()
    tiledlayout(1,3)
    for ii = 1:3
        nexttile
        B_positions_struct.contourf(B(:,ii))
        colorbar()
    end
end
end