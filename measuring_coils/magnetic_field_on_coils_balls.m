home
clear
close all
init_eidors()

addpath('..')
%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 1.5*phantom.radius;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.0025;
phantom.maxsz = 0.005;
phantom.background = 0.503; % muscle at 1 MHz
phantom.extra_format = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
phantom.ball = 0.136; % inflated lung at 1 MHz

current_ampl = 10e-3;
freq = 1e6;

z = linspace(0, phantom.height, 7);
target_centers = zeros(length(z)-2,3);
target_centers(:,3) = z(2:end-1);
target_radius = 0.015;

models = struct('eit', [], 'target_position', [], 'target_radius', []);
tic

phantom.extra = {'', ''};
    
models(1).target_position = [];
models(1).target_radius = [];
models(1).eit = EIT(phantom, current_ampl);
models(1).eit.show_fem()
models(1).eit.calc_elem_current();

for ii = 2:size(target_centers, 1)+1
    phantom.extra = phantom.extra_format(target_radius, target_centers(ii-1,:));
    
    models(ii).target_position = target_centers(ii-1,:);
    models(ii).target_radius = target_radius;
    models(ii).eit = EIT(phantom, current_ampl);
    models(ii).eit.show_fem()
    models(ii).eit.calc_elem_current();
end
toc

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
coil_detectors(1).center = [phantom.radius + 0.01, 0, phantom.elec_vert_position];

%% vertical perpendicolar coils
coil_detectors(2).name = 'Vertical perpendicolar coils';
coil_detectors(2).radius = coil_detectors(1).radius; % m
coil_detectors(2).orientation = [pi/2,0, 0];
coil_detectors(2).center = [phantom.radius + 0.01 + coil_detectors(1).radius, 0, phantom.elec_vert_position];

%% horizontal coplanar coils
coil_detectors(3).name = 'horizontal coplanar coils';
coil_detectors(3).radius = coil_detectors(1).radius; % m
coil_detectors(3).orientation = [0, 0, pi/2];
coil_detectors(3).center = [phantom.radius + 0.01 + coil_detectors(1).radius, 0, phantom.elec_vert_position];

%%
for ii = 1:length(coil_detectors)
    coil_detectors(ii).mother_coil = Coil(coil_detectors(ii).center, ...
                                        coil_detectors(ii).radius, ...
                                        coil_detectors(ii).orientation);
    coil_detectors(ii).coil_system = CoilSystem(coil_detectors(ii).mother_coil, ...
                                        phantom.n_elec);

    % show model with coils
    figure(500+ii)
    clf
    hold on
    coil_detectors(ii).coil_system.show();
    
    show_fem(models(1).eit.img)
    hold off
end

%%
tic
integral_values = zeros(phantom.n_elec*phantom.n_elec, length(coil_detectors), length(models));
for model = 1:length(models)
    for coil = 1:length(coil_detectors)
        values_ = coil_detectors(coil).coil_system.calc_coil_integrals( ...
            models(model).eit, ...
            1:phantom.n_elec);
        integral_values(:, coil, model) = values_(:);
    end
end
magnetic_voltages = 2*pi*freq*integral_values;
toc
%%

% magnetic_voltages = 2*pi*freq*reshape(integral_values, [], 3);
figure
tiledlayout(3,1)

for ii = 1:3
    nexttile
    hold on
    for model = 1:length(models) 
        plot(magnetic_voltages(:,ii, model))
    end
    xlabel('coil index')
    ylabel('coil voltage / V')
    title(coil_detectors(ii).name)
end


% legend('tangential coils', 'perpendicolar coils', 'coplanar coils')

%%

norm_0 = vecnorm(magnetic_voltages(:,:,1), 2, 1);
norms = vecnorm(magnetic_voltages-magnetic_voltages(:,:,1), 2, 1)./norm_0;
norms = squeeze(norms(:,:,2:end));

figure; 
bar(norms)



