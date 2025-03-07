home
clear
close all
init_eidors()

addpath('..')

%%

%%
Dx = [-10,10];
Dy = [-0.5,0.5];
Dz = Dy;

shape_str = [sprintf('solid top    = plane( 0, 0, %f; 0, 0, 1);\n', max(Dz)) ...
          sprintf('solid bot    = plane( 0, 0, %f; 0, 0,-1);\n', min(Dz)) ...
          sprintf('solid xmax   = plane( %f, 0, 0; 1, 0, 0);\n', max(Dx)) ...
          sprintf('solid xmin   = plane( %f, 0, 0;-1, 0, 0);\n', min(Dx)) ...
          sprintf('solid ymax   = plane( 0, %f, 0; 0, 2, 0);\n', max(Dy)) ...
          sprintf('solid ymin   = plane( 0, %f, 0; 0,-2, 0);\n', min(Dy)) ...
          'solid mainobj= top and bot and xmax and xmin and ymax and ymin -maxh=0.5;'];
elec_pos = [min(Dx), 0,  0,   1,  0,  0;
          max(Dx),  0,  0,   -1, 0,  0];
elec_shape=[6,2,0.05];
elec_obj = {'xmin','xmax'};

fwd_model = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

figure(6454)
clf
show_fem(fwd_model, [0,1.012]);
axis equal; view(3);
xlabel('x')
ylabel('y')
zlabel('z')


%%
current_ampl = 1;
stim_pattern = mk_stim_patterns(2, 1, '{op}', '{op}', {'meas_current'}, current_ampl);
            
fwd_model.stimulation = stim_pattern;
[fwd_model.electrode(:).z_contact] = deal(1e-5);


conductivity = 1;
img = mk_image(fwd_model,conductivity);
img.fwd_solve.get_all_meas = 1;

volt_strct = fwd_solve(img);

%%
elem_currents = calc_elem_current(img, volt_strct.volt(:,1));
elem_centers = interp_mesh(fwd_model, 0); % center of elements

elem_volumes = helpers.calc_element_volume(fwd_model.elems, fwd_model.nodes);

%%
B_position = [100, 5, 5];

figure(6454)
clf
show_fem(fwd_model, [0,1.012]);
axis equal; view(3);
hold on
plot3(B_position(1), B_position(2), B_position(3), 'sm')

B_fem = helpers.calc_B_at_points(B_position, elem_centers, elem_currents, elem_volumes)




%% analytical equation
mu0 = 4*pi*1e-7;
r = B_position-[mean(Dx),mean(Dy),mean(Dz)];
dl = current_ampl*[diff(Dx),0,0];
theoretical_values = mu0/(4*pi)*cross(dl,r)./vecnorm(r,2,2).^3

