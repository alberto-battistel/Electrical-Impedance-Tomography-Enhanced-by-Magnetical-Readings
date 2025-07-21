home
clear
close all
%%

init_eidors()

addpath('..')
addpath('progressbar/')

%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 0.02;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001;
phantom.maxsz = 0.002; % 0.005
phantom.background = 1; % muscle at 1 MHz

current_ampl = 10e-3;

target_radius = [0.01];
target_values = [0.1];
target_centers = [0,-0.085,phantom.elec_vert_position];

vert_displacement = 2e-3;
radial_displacement = 1e-3;

% target_model_type = 'right'; % 'simple'
target_model_type = 'right';

%%
blank_model = quick_functions.make_model(phantom, current_ampl);
blank_model = quick_functions.solve_voltage_current(blank_model);

%%
B_positions = (phantom.radius + radial_displacement)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2);
B_positions = [real(B_positions)', imag(B_positions)', phantom.elec_vert_position*ones(phantom.n_elec,1)+ vert_displacement];


%%
B0 = quick_functions.calc_B(blank_model, B_positions);


%%

all_target_values = combinations(target_radius, target_values);

target_models = cell(size(all_target_values, 1), 1);


for ii = 1:size(all_target_values, 1)
    target_models{ii} = quick_functions.mk_model_target(blank_model, ...
        target_centers, ...
        all_target_values.target_radius(ii), ...
        all_target_values.target_values(ii), target_model_type);
end

%%
figure(1)
clf
hold on
show_fem(target_models{1}.img)
children = gca().Children;
children(end).EdgeAlpha = 0.4;
set(gca, "CameraPosition", [0 -1.3 0.5])

plot3(B_positions(:,1),B_positions(:,2),B_positions(:,3), 'om')

for ii = 1:4
    text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii))
end

xlabel('x / m')
ylabel('y / m')
zlabel('z / m')

%%
progressbar('Solve Target Models')
for mdl = 1:length(target_models)
    target_models{mdl} = quick_functions.solve_voltage_current(target_models{mdl});
    progressbar(mdl/length(target_models))
end

%%
B_targets = zeros(size(B0,1), size(B0,2), length(target_models));

progressbar(0,0)
for mdl = 1:length(target_models)
    B_targets(:,:,:,mdl) = quick_functions.calc_B(target_models{mdl}, B_positions);
    progressbar(mdl/length(target_models), []) % 
end




%%
values = B_targets;
values_0 = B0;
fun = @(y, y0) vecnorm(y-y0,2,1)./vecnorm(y0,2,1);

normalize_diff = zeros(length(target_models),3);

figure(100)
tiledlayout(2,3)
for ii = 1:3
    nexttile
    hold on
    plot(squeeze(values(:,ii,:)))
    plot(values_0(:,ii), 'k')
    hold off

    normalize_diff(:,ii) = fun(squeeze(values(:,ii,mdl)), values_0(:,ii));
    subtitle(sprintf("norm diff = %.3g", normalize_diff(:,ii)))
end

for ii = 1:3
    nexttile
    plot(squeeze(values(:,ii,:)) - values_0(:,ii))

end

% it may do strange things if you have more than 1 target
normalize_diff_EIT = fun(target_models{1}.volt_strct.meas, blank_model.volt_strct.meas);
values_rms_EIT = rms(target_models{1}.volt_strct.meas);

figure(500)
tiledlayout(2,1)
nexttile
hold on
plot(target_models{1}.volt_strct.meas, 'r')
plot(blank_model.volt_strct.meas, 'k')
hold off
ylabel('V / V')
subtitle(sprintf("norm diff = %.3g", normalize_diff_EIT))

nexttile
plot(target_models{1}.volt_strct.meas-blank_model.volt_strct.meas, 'k')
xlabel('Measurement Index')
ylabel('\DeltaV / V')


















