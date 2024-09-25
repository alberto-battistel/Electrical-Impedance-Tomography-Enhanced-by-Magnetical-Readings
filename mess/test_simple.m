% taken from https://eidors3d.sourceforge.net/doc/index.html?eidors/models/interp_mesh.html
% and https://eidors3d.sourceforge.net/doc/index.html?eidors/graphics/matlab/show_current.html
% do_unit_test

close all
clear

mu0 = 4*pi*1e-7;

fmdl.nodes= [4,6,8,4,6,8;2,2,2,5,5,5]';
fmdl.elems=[1,2,4;2,4,5;2,3,5;3,5,6];
fmdl.type='fwd_model';
elem_centers = interp_mesh(fmdl);
show_fem(fmdl); 
hold on ;
for i=1:size(elem_centers,3)
    plot(elem_centers(:,1,i), elem_centers(:,2,i),'*'); 
end
hold off


fmdl.electrode(1).nodes = [1,4]; fmdl.electrode(1).z_contact = 0.01;
fmdl.electrode(2).nodes = [3,6]; fmdl.electrode(2).z_contact = 0.01;
fmdl.gnd_node = 1;
fmdl.stimulation(1).stim_pattern = [1e-3;-1e-3]; % 1 mA
fmdl.stimulation(1).meas_pattern = [1,-1];
fmdl.solve = @fwd_solve_1st_order;
fmdl.system_mat = @system_mat_1st_order;
fmdl.type = 'fwd_model';
fmdl.normalize_measurements= 0;
img = mk_image(fmdl,1); 
img.fwd_solve.get_all_meas = 1;

% show_current(img); % it is as e_curr (e_curr = show_current(img))

vh = fwd_solve(img);
e_curr = calc_elem_current(img, vh.volt(:,1));

elem_areas = calc_area(fmdl.nodes, fmdl.elems);

point = [1,1];

r = point - elem_centers;

r_mag = norm(r);

dB = mu0/(4*pi)*cross_prod_2d(e_curr.*elem_areas, r)./norm(r.^3); % magnetic flux density, in T

%%
function axb = cross_prod_2d(a,b)
    axb = a(:,1).*b(:,2) - a(:,2).*b(:,1);
end