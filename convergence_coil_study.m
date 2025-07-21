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
phantom.radius = 0.05; % 0.1
phantom.height = 0.04; % 0.04
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001; % 0.001
phantom.maxsz = 0.005; % 0.005
phantom.background = 1; % muscle at 1 MHz

current_ampl = 10e-3;

maxsz = [0.01, 0.008, 0.0075, 0.0025];


magnetic_sensor_info.radial_displacement = 10e-3; % mm
magnetic_sensor_info.vert_displacement = 10e-3; % mm
magnetic_sensor_info.ang_displacement = 0; % degree
magnetic_sensor_info.coil_radius = 5e-3;

folder_name = "convergence_study_"+string(datetime('now','format', 'yyyyMMdd_HH_mm_ss'));

%%
models = cell(length(maxsz),2);

progressbar('Make models')
for mdl = 1:length(maxsz)
    phantom.maxsz = maxsz(mdl);
    models{mdl,1} = quick_functions.make_model(phantom, current_ampl);
    models{mdl,2} = models{mdl,1};
    models{mdl,2}.img.fwd_model.solve = @fwd_solve_higher_order;
    models{mdl,2}.img.fwd_model.system_mat = @system_mat_higher_order;
    models{mdl,2}.img.fwd_model.approx_type = 'tet10'; %Quadratic
    progressbar(mdl/length(maxsz))
end

%%
target_models = cell(length(maxsz),2);
target_centers = [0, -phantom.radius/2, phantom.elec_vert_position];
target_radius = 0.01;
target_value = 0.1;
target_model_type = 'right';

progressbar('Make target models')
for mdl = 1:length(maxsz)
    target_models{mdl,1} = quick_functions.mk_model_target(models{mdl,1}, ...
        target_centers, ...
        target_radius, ...
        target_value, ...
        target_model_type);
    target_models{mdl,2} = target_models{mdl,1};
    target_models{mdl,2}.img.fwd_model.solve = @fwd_solve_higher_order;
    target_models{mdl,2}.img.fwd_model.system_mat = @system_mat_higher_order;
    target_models{mdl,2}.img.fwd_model.approx_type = 'tet10'; %Quadratic
    progressbar(mdl/length(maxsz))
end

%%
empty_models = cell(length(maxsz),2);
target_value = phantom.background;

progressbar('Make empty target models')
for mdl = 1:length(maxsz)
    empty_models{mdl,1} = quick_functions.mk_model_target(models{mdl,1}, ...
        target_centers, ...
        target_radius, ...
        target_value, ...
        target_model_type);
    empty_models{mdl,2} = empty_models{mdl,1};
    empty_models{mdl,2}.img.fwd_model.solve = @fwd_solve_higher_order;
    empty_models{mdl,2}.img.fwd_model.system_mat = @system_mat_higher_order;
    empty_models{mdl,2}.img.fwd_model.approx_type = 'tet10'; %Quadratic
    progressbar(mdl/length(maxsz))
end

%%
progressbar('Solve models')
for mdl = 1:length(maxsz)
    models{mdl,1} = quick_functions.solve_voltage_current(models{mdl,1});
    models{mdl,2} = quick_functions.solve_voltage_current(models{mdl,2});
    progressbar(mdl/length(maxsz))
end

progressbar('Solve target models')
for mdl = 1:length(maxsz)
    target_models{mdl,1} = quick_functions.solve_voltage_current(target_models{mdl,1});
    target_models{mdl,2} = quick_functions.solve_voltage_current(target_models{mdl,2});
    progressbar(mdl/length(maxsz))
end

progressbar('Solve empty target models')
for mdl = 1:length(maxsz)
    empty_models{mdl,1} = quick_functions.solve_voltage_current(empty_models{mdl,1});
    empty_models{mdl,2} = quick_functions.solve_voltage_current(empty_models{mdl,2});
    progressbar(mdl/length(maxsz))
end

%%
V0 = cell(size(models));
for mdl = 1:length(maxsz)
    V0{mdl,1} = models{mdl,1}.volt_strct.meas;
    V0{mdl,2} = models{mdl,2}.volt_strct.meas;
end

V = cell(size(models));
for mdl = 1:length(maxsz)
    V{mdl,1} = target_models{mdl,1}.volt_strct.meas;
    V{mdl,2} = target_models{mdl,2}.volt_strct.meas;
end

Ve = cell(size(models));
for mdl = 1:length(maxsz)
    Ve{mdl,1} = empty_models{mdl,1}.volt_strct.meas;
    Ve{mdl,2} = empty_models{mdl,2}.volt_strct.meas;
end

%%
B_pos = (phantom.radius + magnetic_sensor_info.radial_displacement)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2 +1j*magnetic_sensor_info.ang_displacement/180*pi);
B_positions = [real(B_pos)', imag(B_pos)', phantom.elec_vert_position*ones(phantom.n_elec,1)+ magnetic_sensor_info.vert_displacement];

coil_info = struct("radius", magnetic_sensor_info.coil_radius, ...
                    "coil_minimal_distance", 0, ...
                    "displacement", [magnetic_sensor_info.radial_displacement, 0, magnetic_sensor_info.vert_displacement]);

coils = quick_functions.make_coil_systems(coil_info, phantom);


%%

long_coil_0 = cell(size(models));
progressbar('Calculate coil 0')
for mdl = 1:length(maxsz)
    for ii = 1:3
        integral_values = coils{ii,1}.calc_coil_integrals(models{mdl,1}, 1:16);
        long_coil_0{mdl,1}(:,ii) = integral_values(:);
        integral_values = coils{ii,1}.calc_coil_integrals(models{mdl,2}, 1:16);
        long_coil_0{mdl,2}(:,ii) = integral_values(:);
    end
    progressbar(mdl/length(maxsz))
end

long_coil_t = cell(size(target_models));
progressbar('Calculate coil target')
for mdl = 1:length(maxsz)
    for ii = 1:3
        integral_values = coils{ii,1}.calc_coil_integrals(target_models{mdl,1}, 1:16);
        long_coil_t{mdl,1}(:,ii) = integral_values(:);
        integral_values = coils{ii,1}.calc_coil_integrals(target_models{mdl,2}, 1:16);
        long_coil_t{mdl,2}(:,ii) = integral_values(:);
    end
    progressbar(mdl/length(maxsz))
end

long_coil_e = cell(size(target_models));
progressbar('Calculate coil empty target')
for mdl = 1:length(maxsz)
    for ii = 1:3
        integral_values = coils{ii,1}.calc_coil_integrals(empty_models{mdl,1}, 1:16);
        long_coil_e{mdl,1}(:,ii) = integral_values(:);
        integral_values = coils{ii,1}.calc_coil_integrals(empty_models{mdl,2}, 1:16);
        long_coil_e{mdl,2}(:,ii) = integral_values(:);
    end
    progressbar(mdl/length(maxsz))
end


%%
figure(1)
tiledlayout(1,3)
nexttile
hold on
show_fem(models{2,2}.img)

children = gca().Children;
children(end).EdgeAlpha = 0.4;

set(gca, "CameraPosition", [0 -1.3 0.5])

plot3(B_positions(:,1),B_positions(:,2),B_positions(:,3), 'om')

for ii = [1:9]
    text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii), 'color', 'r')
end
xlabel('x / m')
ylabel('y / m')
zlabel('z / m')

nexttile
hold on
show_fem(target_models{2,2}.img)

children = gca().Children;
children(end).EdgeAlpha = 0.4;

set(gca, "CameraPosition", [0 -1.3 0.5])

plot3(B_positions(:,1),B_positions(:,2),B_positions(:,3), 'om')

for ii = [1:9]
    text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii), 'color', 'r')
end
xlabel('x / m')
ylabel('y / m')
zlabel('z / m')

nexttile
hold on
show_fem(empty_models{2,2}.img)

children = gca().Children;
children(end).EdgeAlpha = 0.4;

set(gca, "CameraPosition", [0 -1.3 0.5])

plot3(B_positions(:,1),B_positions(:,2),B_positions(:,3), 'om')

for ii = [1:9]
    text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii), 'color', 'r')
end
xlabel('x / m')
ylabel('y / m')
zlabel('z / m')

%%
k0 = [1:17:256, 2:17:256, 0:17:256];
k0 = sort(k0(k0>0));
idx = 1:256; 
ll = ones(256,1); 
ll(k0) = 0;
ll=logical(ll);

for mdl = 1:size(long_coil_0,1)
    for ii = 1:size(V0,2)
        coil_0{mdl,ii} = long_coil_0{mdl,ii}(ll,:);
        coil_t{mdl,ii} = long_coil_t{mdl,ii}(ll,:);
        coil_e{mdl,ii} = long_coil_e{mdl,ii}(ll,:);
    end
end

figure(2)
tiledlayout(3,1)
for ii = 1:3
    nexttile
    hold on
    plot(idx, long_coil_0{end,end}(:,ii))
    plot(idx(ll), coil_0{end,end}(:,ii), '.')
end

%%

figure(10)
tiledlayout(1,2)
for mdl = 1:length(maxsz)
    nexttile(1)
    hold on
    plot(V0{mdl,1})
    nexttile(2)
    hold on
    plot(V0{mdl,2})
end

figure(11)
tiledlayout(1,2)
for mdl = 1:length(maxsz)
    nexttile(1)
    hold on
    plot(V{mdl,1})
    nexttile(2)
    hold on
    plot(V{mdl,2})
end

figure(12)
tiledlayout(1,2)
for mdl = 1:length(maxsz)
    nexttile(1)
    hold on
    plot(V{mdl,1}-V0{mdl,1})
    nexttile(2)
    hold on
    plot(V{mdl,2}-V0{mdl,2})
end

figure(13)
tiledlayout(1,2)
for mdl = 1:length(maxsz)
    nexttile(1)
    hold on
    plot(Ve{mdl,1}-V0{mdl,1})
    nexttile(2)
    hold on
    plot(Ve{mdl,2}-V0{mdl,2})
end

%%
figure(20)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(coil_0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(coil_0{mdl,2}(:,ii))
    end
end

figure(21)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(coil_t{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(coil_t{mdl,2}(:,ii))
    end
end


figure(22)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(coil_t{mdl,1}(:,ii)-coil_0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(coil_t{mdl,2}(:,ii)-coil_0{mdl,2}(:,ii))
    end
end

figure(23)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(coil_e{mdl,1}(:,ii)-coil_0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(coil_e{mdl,2}(:,ii)-coil_0{mdl,2}(:,ii))
    end
end

%%
reference = V0{end,2};
norm_diff_V = zeros(size(V0));
for mdl = 1:size(V0,1)
    for ll = 1:size(V0,2)
        norm_diff_V(mdl,ll) = vecnorm(V0{mdl,ll}-reference,2,1)./vecnorm(reference,2,1);
    end
end

figure(100)
plot(maxsz, norm_diff_V)



%%
reference = coil_0{end,2};
norm_diff_B = zeros(size(V0,1), size(V0,2), 3);
for mdl = 1:size(coil_0,1)
    for ll = 1:size(V0,2)
        norm_diff_B(mdl,ll,:) = vecnorm(coil_0{mdl,ll}-reference,2,1)./vecnorm(reference,2,1);
    end
end

figure(200)
tiledlayout(3,1);
for ii = 1:3
    nexttile
    plot(maxsz, norm_diff_B(:,:,ii))
end

%%
norm_diff_VV0 = zeros(size(V0));
norm_diff_V0Ve = zeros(size(V0));
for mdl = 1:size(V0,1)
    for ll = 1:size(V0,2)
        norm_diff_VV0(mdl,ll) = vecnorm(V0{mdl,ll}-V{mdl,ll},2,1)./vecnorm(V0{mdl,ll},2,1);
        norm_diff_V0Ve(mdl,ll) = vecnorm(V0{mdl,ll}-Ve{mdl,ll},2,1)./vecnorm(V0{mdl,ll},2,1);
    end
end

figure(300)
tiledlayout(1,2)
nexttile
plot(maxsz, norm_diff_VV0)
nexttile
plot(maxsz, norm_diff_V0Ve)

%%
norm_diff_BB0 = zeros(size(V0,1), size(V0,2), 3);
norm_diff_B0Be = zeros(size(V0,1), size(V0,2), 3);

for mdl = 1:size(coil_0,1)
    for ll = 1:size(coil_0,2)
        norm_diff_BB0(mdl,ll,:) = vecnorm(coil_0{mdl,ll}-coil_t{mdl,ll},2,1)./vecnorm(coil_0{mdl,ll},2,1);
        norm_diff_B0Be(mdl,ll,:) = vecnorm(coil_0{mdl,ll}-coil_e{mdl,ll},2,1)./vecnorm(coil_0{mdl,ll},2,1);
    end
end

figure(301)
t= tiledlayout(3,2);
for ii = 1:3
    nexttile(tilenum(t,ii,1))
    plot(maxsz, norm_diff_BB0(:,:,ii))
end
for ii = 1:3
    nexttile(tilenum(t,ii,2))
    plot(maxsz, norm_diff_B0Be(:,:,ii))
end
%%

save_all_figures(folder_name)

if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end

cwd = pwd;
cd(folder_name)
save("data.mat", ...
    "long_coil_0", ...
    "long_coil_t", ...
    "long_coil_e", ...
    "coil_0", ...
    "coil_t", ...
    "coil_e", ...
    "coils", ...
    "idx", ...
    "V0", ...
    "maxsz", ...
    "V", ...
    "Ve", ...
    "phantom", ...
    "magnetic_sensor_info",...
    "norm_diff_V", ...
    "norm_diff_B", ...
    "norm_diff_VV0", ...
    "norm_diff_V0Ve", ...
    "norm_diff_BB0", ...
    "norm_diff_B0Be")
cd(cwd)

