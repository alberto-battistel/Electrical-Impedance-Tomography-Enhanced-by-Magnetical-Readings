function coil_detectors = make_coil_detectors(phantom_radius, phantom_elec_vert_position, n_coils)
%%
coil_detectors = struct('name', [], ...
                        'center', [], ...
                        'radius', [], ...
                        'orientation', [], ...
                        'mother_coil', [], ...
                        'coil_system', []);

%% horizontal tangential coils
coil_detectors(1).name = 'Horizontal tangential coils';
coil_detectors(1).radius = 0.005; % m
coil_detectors(1).orientation = [0,pi/2, 0];
coil_detectors(1).center = [phantom_radius + 0.01, 0, phantom_elec_vert_position];

%% vertical perpendicolar coils
coil_detectors(2).name = 'Vertical perpendicolar coils';
coil_detectors(2).radius = coil_detectors(1).radius; % m
coil_detectors(2).orientation = [pi/2,0, 0];
coil_detectors(2).center = [phantom_radius + 0.01 + coil_detectors(1).radius, 0, phantom_elec_vert_position];

%% horizontal coplanar coils
coil_detectors(3).name = 'horizontal coplanar coils';
coil_detectors(3).radius = coil_detectors(1).radius; % m
coil_detectors(3).orientation = [0, 0, pi/2];
coil_detectors(3).center = [phantom_radius + 0.01 + coil_detectors(1).radius, 0, phantom_elec_vert_position];


for ii = 1:length(coil_detectors)
    coil_detectors(ii).mother_coil = Coil(coil_detectors(ii).center, ...
                                        coil_detectors(ii).radius, ...
                                        coil_detectors(ii).orientation);
    coil_detectors(ii).coil_system = CoilSystem(coil_detectors(ii).mother_coil, ...
                                        n_coils);
end

end