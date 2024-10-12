% use eidors to make the mesh of a disk
imdl = mk_common_model('a2c2',8);
fmdl = imdl.fwd_model;

% crop
fmdl = crop_model(fmdl,inline('x<0.5','x','y','z'));
fmdl = crop_model(fmdl,inline('y<-0.2','x','y','z'));
figure(1)
show_fem(fmdl)

original_points = [fmdl.nodes, zeros(length(fmdl.nodes),1)];
connectivity_list = fmdl.elems;


[original_areas] = helpers.calc_area(original_points, connectivity_list);
function_values = [1:length(original_points);...
                    10+(1:length(original_points));...
                    20+(1:length(original_points))]';
% function_values = [1:length(original_points)]';

mesh_with_values = struct('points', original_points, 'connectivity_list', connectivity_list, 'areas', original_areas,'values', function_values);

final_integral = integral_on_mesh(mesh_with_values)

