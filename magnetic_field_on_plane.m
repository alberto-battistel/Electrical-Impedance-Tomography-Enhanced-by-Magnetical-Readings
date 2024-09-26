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

el_pos = [-360/n_elec/2+(0:n_elec-1).'/n_elec*360,elec_vert_position.*ones(16,1)];
el_sz  = [0.01,0,0.05].*ones(size(el_pos,1),3);

fmdl = ng_mk_cyl_models([phantom_height,phantom_radius,maxsz], el_pos, el_sz);
imdl = mk_common_model('a2c2',16); % Will replace most 
imdl.fwd_model = fmdl;
imdl.normalize_measurements = 0;
stim_pattern = mk_stim_patterns(size(el_pos,1),1,'{ad}','{ad}',{},10e-3);

imdl.fwd_model.stimulation = stim_pattern;

img_h = mk_image(imdl, 0.503); % muscle cond 1 MHz
img_h.fwd_solve.get_all_meas = 1;

figure()
show_fem(fmdl,[0,1.012])

%%
%
elem_centers = interp_mesh(imdl.fwd_model, 0); % center of elements
elem_areas = helpers.calc_area(imdl.fwd_model.nodes, imdl.fwd_model.elems);

vh = fwd_solve(img_h);

figure()
show_current(img_h,vh.volt(:,1));

e_curr = calc_elem_current(img_h, vh.volt(:,1));

%% vert plane tangential to phantom
n_points = 50;
x_dir = linspace(-phantom_radius, phantom_radius, n_points);
y_dir = phantom_radius;
z_dir = linspace(0, phantom_height, n_points);
[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions = [xx(:), yy(:), zz(:)];

figure()
clf
hold on
show_fem(fmdl,[0,1.012])
plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), 'ro','MarkerFaceColor','r')

B = helpers.calc_B_at_points(B_positions, elem_centers, e_curr, elem_areas);

%%
abs_B = vecnorm(B,2,2);
[xx, zz] = meshgrid(x_dir, z_dir);

figure()
surf(xx,zz,reshape(abs_B,size(xx)))

%%
figure()
tiledlayout(1,3)
for ii = 1:3
    nexttile
    contourf(xx, zz, reshape(B(:,ii),size(xx)))
end


%% horz plane tangential to phantom
n_points = 50;
x_dir = linspace(-phantom_radius, phantom_radius, n_points);
y_dir = linspace(phantom_radius, 3*phantom_radius, n_points);
z_dir = phantom_height/2;
[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions = [xx(:), yy(:), zz(:)];

figure()
clf
hold on
show_fem(fmdl,[0,1.012])
plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), 'ro','MarkerFaceColor','r')

B = helpers.calc_B_at_points(B_positions, elem_centers, e_curr, elem_areas);

%%
abs_B = vecnorm(B,2,2);

[xx, zz] = meshgrid(x_dir, y_dir);

figure()
surf(xx,zz,reshape(abs_B,size(xx)))

%%
figure()
tiledlayout(1,3)
for ii = 1:3
    nexttile
    contourf(xx, zz, reshape(B(:,ii),size(xx)))
end


%% vert plane sagital to phantom
n_points = 50;
x_dir = 0;
y_dir = linspace(phantom_radius, 3*phantom_radius, n_points);
z_dir = linspace(0, phantom_height, n_points);
[xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
B_positions = [xx(:), yy(:), zz(:)];

figure()
clf
hold on
show_fem(fmdl,[0,1.012])
plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), 'ro','MarkerFaceColor','r')

B = helpers.calc_B_at_points(B_positions, elem_centers, e_curr, elem_areas);

%%
abs_B = vecnorm(B,2,2);

[xx, zz] = meshgrid(z_dir, y_dir);

figure()
surf(xx,zz,reshape(abs_B,size(xx)))

%%
figure()
tiledlayout(1,3)
for ii = 1:3
    nexttile
    contourf(xx, zz, reshape(B(:,ii),size(xx)))
end