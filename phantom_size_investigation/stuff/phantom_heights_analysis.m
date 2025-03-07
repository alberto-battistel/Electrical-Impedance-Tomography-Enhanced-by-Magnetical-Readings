home
clear
close all
init_eidors()

addpath('..')
%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = [];
phantom.elec_vert_position = [];
phantom.max_el_sz = 0.001; %0.0025
phantom.maxsz = 0.0025;
phantom.background = 0.503; % muscle at 1 MHz
phantom.extra = {'', ''};
phantom.extra_format = {'', ''};

n_coils = 1;

current_ampl = 10e-3;
freq = 1e6;

phantom_heights = 0.1:0.1:0.5;


models = struct('eit', [], 'phantom_height', [], 'coil_detectors', []);

tic
for ii = 1:length(phantom_heights)
    phantom.height = phantom_heights(ii);
    phantom.elec_vert_position = phantom.height/2;

    fprintf('Phantom height: %.3f m\n', phantom.height)
    
    models(ii).phantom_height = phantom.height;

    models(ii).eit = EIT(phantom, current_ampl);
end


%%

for model = 1:length(models)
    models(model).eit.calc_elem_current(1);
end


%%
for ii = 1:length(phantom_heights)
    models(ii).coil_detectors = make_coil_detectors(phantom.radius, ...
        phantom_heights(ii)/2, n_coils);
end



%% plot 
% strange

% for model = 1:length(models)
%     for iii = 1:length(models(model).coil_detectors)
%         % show model with coils
%         figure(100*model+iii)
%         clf
%         hold on
%         models(model).coil_detectors(iii).coil_system.show();
% 
%         % show_fem(models(ii).eit.img)
%         hold off
%     end
% end

%%

integral_values = zeros(n_coils, length(models(1).coil_detectors), length(models));

for model = 1:length(models)
    for coil = 1:length(models(model).coil_detectors)
        values_ = models(model).coil_detectors(coil).coil_system.calc_coil_integrals( ...
            models(model).eit, ...
            1);
        integral_values(:, coil, model) = values_(:);
    end
end


%%
magnetic_voltages = squeeze(2*pi*freq*integral_values);

norm_magnetic_voltages = squeeze(magnetic_voltages./magnetic_voltages(:,end));
abs_rel_diff_magnetic_voltages = abs(norm_magnetic_voltages(:,1:end-1)-1);

%%

figure(100)
tiledlayout(3,1)

x = [models(1:end-1).phantom_height];
for comp = 1:3
    nexttile
    semilogy(x, abs_rel_diff_magnetic_voltages(comp,:))
end


%%

figure(200)
tiledlayout(3,1)

x = [models(:).phantom_height];
for comp = 1:3
    nexttile
    plot(x, magnetic_voltages(comp,:))
end

toc
%%


