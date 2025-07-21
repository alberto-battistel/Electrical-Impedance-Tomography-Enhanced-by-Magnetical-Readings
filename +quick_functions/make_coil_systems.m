function coils = make_coil_systems(coil_info, phantom)

coils = cell(3,1);

%% tangential coils
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance, 0, phantom.elec_vert_position];
coil_info.center = coil_info.center + coil_info.displacement;
coil_info.radius = coil_info.radius;
coil_info.orientation = [0, pi/2, 0];

coils{1} = quick_functions.make_coil(coil_info);

%% axial coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.center = coil_info.center + coil_info.displacement;
coil_info.orientation = [pi/2,0, 0];

coils{2} = quick_functions.make_coil(coil_info);

%% coplanar coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.center = coil_info.center + coil_info.displacement;
coil_info.orientation = [0, 0, pi/2];

coils{3} = quick_functions.make_coil(coil_info);

end