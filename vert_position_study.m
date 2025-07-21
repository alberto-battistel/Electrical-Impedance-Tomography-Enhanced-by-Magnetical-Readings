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
phantom.height = 0.04;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001; % 0.001
phantom.maxsz = 0.0025; % 0.005
phantom.background = 1; % muscle at 1 MHz

current_ampl = 10e-3;

sensor_vert_position = [1e-3, 2e-3, 3e-3, 5e-3];

target_radius = [0.01];
target_values = [0.1];
target_centers = [0,-0.085,phantom.elec_vert_position];

all_target_values = combinations(target_radius, target_values);

% target_model_type = 'right'; % 'simple'
target_model_type = 'right';
folder_name = "vert_position_study_"+string(datetime('now','format', 'yyyyMMdd_HH_mm_ss'));

%%
blank_model = quick_functions.make_model(phantom, current_ampl);
blank_model = quick_functions.solve_voltage_current(blank_model);

%%
target_models = cell(size(all_target_values, 1), 1);

progressbar('Make Target Models')
for ii = 1:size(all_target_values, 1)
    target_models{ii} = quick_functions.mk_model_target(blank_model, ...
        target_centers, ...
        all_target_values.target_radius(ii), ...
        all_target_values.target_values(ii), ...
        target_model_type);
    progressbar(mdl/length(target_models))
end

%%
progressbar('Solve Target Models')
for mdl = 1:length(target_models)
    target_models{mdl} = quick_functions.solve_voltage_current(target_models{mdl});
    progressbar(mdl/length(target_models))
end

%%
B_positions = cell(length(sensor_vert_position),1);
B0 = cell(length(sensor_vert_position),1);
B_targets = cell(length(sensor_vert_position),1);


for pos = 1:length(sensor_vert_position)

    vert_displacement = sensor_vert_position(pos);
    B_pos = (phantom.radius + 0.001)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2);
    B_positions{pos,1} = [real(B_pos)', imag(B_pos)', phantom.elec_vert_position*ones(phantom.n_elec,1)+ vert_displacement];
        
    B0{pos,1} = quick_functions.calc_B(blank_model, B_positions{pos,1});

    B_targets{pos,1} = zeros(size(B0{pos,1},1), size(B0{pos,1},2), length(target_models));

    progressbar(sprintf('Calc B Target Models pos %d', pos))
    for mdl = 1:length(target_models)
        B_targets{pos,1}(:,:,:,mdl) = quick_functions.calc_B(target_models{mdl}, B_positions{pos,1});
        progressbar(mdl/length(target_models)) % 
    end
end

%%
figure(1)
tiledlayout(2,1)
nexttile
hold on
show_fem(target_models{1}.img)

children = gca().Children;
children(end).EdgeAlpha = 0.4;

set(gca, "CameraPosition", [0 -1.3 0.5])

plot3(B_positions{1,1}(:,1),B_positions{1,1}(:,2),B_positions{1,1}(:,3), 'om')

for ii = [1:9]
    text(B_positions{1,1}(ii,1), B_positions{1,1}(ii,2), B_positions{1,1}(ii,3), sprintf('%d', ii), 'color', 'r')
end
xlabel('x / m')
ylabel('y / m')
zlabel('z / m')

nexttile
hold on
show_fem(target_models{1}.img)

children = gca().Children;
children(end).EdgeAlpha = 0.4;

set(gca, "CameraPosition", [0 -1.4 0])
set(gca, "View", [0 0])

plot3(B_positions{1,1}(:,1),B_positions{1,1}(:,2),B_positions{1,1}(:,3), 'om')

ll = 9;
for ii = 1:length(sensor_vert_position)
    plot3(B_positions{ii,1}(ll,1),B_positions{ii,1}(ll,2),B_positions{ii,1}(ll,3), 'om')
end

for pos = 1:length(sensor_vert_position)
    text(B_positions{pos,1}(ll,1)+0.001, B_positions{pos,1}(ll,2), B_positions{pos,1}(ll,3), sprintf('%.0f mm', sensor_vert_position(pos)*1000), 'color', 'r')
end

xlim([-0.02,0.02])
zlim(phantom.elec_vert_position+[-0.005,0.005])
xlabel('x / m')
ylabel('y / m')
zlabel('z / m')


%%
fun = @(y, y0) vecnorm(y-y0,2,1)./vecnorm(y0,2,1);
normalize_diff_list = cell(length(sensor_vert_position),1);
values_rms_list = cell(length(sensor_vert_position),1);

B_labels = {'\rho', '\theta', 'z'};

for pos = 1:length(sensor_vert_position)

    values = B_targets{pos,1};
    values_0 = B0{pos,1};
   
    
    normalize_diff_list{pos,1} = zeros(length(target_models),3);
    values_rms_list{pos,1} = zeros(length(target_models),3);

    figure(100+pos)
    t1 = tiledlayout(2,3);
    for ii = 1:3
        nexttile
        hold on
        plot(squeeze(values(:,ii,:)), 'r')
        plot(values_0(:,ii), 'k')
        hold off
        
        values_rms_list{pos,1}(1,ii) = rms(squeeze(values(:,ii,:)));

        normalize_diff_list{pos,1}(:,ii) = fun(squeeze(values(:,ii,mdl)), values_0(:,ii));
        title(B_labels{ii})
        % title(sprintf("RMS: %.2g", values_rms_list{pos,1}(ii)))
        % subtitle(sprintf("Diff. = %.3g", normalize_diff_list{pos,1}(:,ii)))
        % ax = gca;
        % ax.TitleHorizontalAlignment = 'right';
        % xlabel('Sensor Index')
        ylabel('B / T')
    end
    title(t1,sprintf('Vert. Pos. %.4f', sensor_vert_position(pos)))

    for ii = 1:3
        nexttile
        plot(squeeze(values(:,ii,:)) - values_0(:,ii), 'k')
        
        % title(B_labels{ii})
        xlabel('Measurement Index')
        ylabel('\Delta B / T')
    end
    

end

% it may do strange things if you have more than 1 target
figure(500)
tiledlayout(2,1)
nexttile
hold on
plot(target_models{1}.volt_strct.meas, 'r')
plot(blank_model.volt_strct.meas, 'k')
hold off
ylabel('V / V')

nexttile
plot(target_models{1}.volt_strct.meas-blank_model.volt_strct.meas, 'k')
xlabel('Measurement Index')
ylabel('\DeltaV / V')

normalize_diff_EIT = fun(target_models{1}.volt_strct.meas, blank_model.volt_strct.meas);
values_rms_EIT = rms(target_models{1}.volt_strct.meas);

%%
% it may do strange things if you have more than 1 target
x = categorical(sensor_vert_position);
normalize_diff = cell2mat(normalize_diff_list);
values_rms = cell2mat(values_rms_list);

figure(1000)
t = tiledlayout(2,1);
nexttile
hold on
bar(x, normalize_diff)
% set(gca, 'yscale', 'log')
plot(x, normalize_diff_EIT.*ones(1,length(sensor_vert_position)), 'dk', 'MarkerFaceColor','k')
xlabel('Vert. Position / m')
ylabel('Norm. Diff.')
legend('\rho', '\theta', 'z', 'EIT')
% xlim([x(1), x(end)])
hold off

nexttile
hold on
bar(x, values_rms)
set(gca, 'yscale', 'log')
% plot(x, values_rms_EIT.*ones(1,length(sensor_vert_position)), 'dk', 'MarkerFaceColor','k')
xlabel('Vert. Position / m')
ylabel('RMS / T')
% xlim([x(1), x(end)])
ylim([1e-10,1e-7])
hold off


%% 
save_all_figures(folder_name)

if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end

cwd = pwd;
cd(folder_name)
save("data.mat", "x", ...
    "normalize_diff", ...
    "values_rms", ...
    "normalize_diff_EIT", ...
    "B0", ...
    "B_targets", ...
    "target_models", ...
    "blank_model")
cd(cwd)
