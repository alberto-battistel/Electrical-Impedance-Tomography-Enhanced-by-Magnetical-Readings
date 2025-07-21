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

maxsz = [0.008, 0.006, 0.004, 0.002];

radial_displacement = 10e-3; % mm
vert_displacement = 10e-3; % mm
ang_displacement = 0; % degree

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
B_pos = (phantom.radius + radial_displacement)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2 +1j*ang_displacement/180*pi);
B_positions = [real(B_pos)', imag(B_pos)', phantom.elec_vert_position*ones(phantom.n_elec,1)+ vert_displacement];

%%
lB0 = cell(size(models));
progressbar('Calculate B0')
for mdl = 1:length(maxsz)
    lB0{mdl,1} = quick_functions.calc_B(models{mdl,1}, B_positions);
    lB0{mdl,2} = quick_functions.calc_B(models{mdl,2}, B_positions);
    progressbar(mdl/length(maxsz))
end

lB = cell(size(target_models));
progressbar('Calculate B')
for mdl = 1:length(maxsz)
    lB{mdl,1} = quick_functions.calc_B(target_models{mdl,1}, B_positions);
    lB{mdl,2} = quick_functions.calc_B(target_models{mdl,2}, B_positions);
    progressbar(mdl/length(maxsz))
end

lBe = cell(size(target_models));
progressbar('Calculate B')
for mdl = 1:length(maxsz)
    lBe{mdl,1} = quick_functions.calc_B(empty_models{mdl,1}, B_positions);
    lBe{mdl,2} = quick_functions.calc_B(empty_models{mdl,2}, B_positions);
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
% ll(k0) = 0;
ll=logical(ll);

for mdl = 1:size(lB0,1)
    for ii = 1:size(V0,2)
        B0{mdl,ii} = lB0{mdl,ii}(ll,:);
        B{mdl,ii} = lB{mdl,ii}(ll,:);
        Be{mdl,ii} = lBe{mdl,ii}(ll,:);
    end
end

figure(2)
clf
hold on
plot(idx, lB0{end,end}(:,end))
plot(idx(ll), B0{end,end}(:,end), '.')

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
        plot(B0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(B0{mdl,2}(:,ii))
    end
end

figure(21)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(B{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(B{mdl,2}(:,ii))
    end
end


figure(22)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(B{mdl,1}(:,ii)-B0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(B{mdl,2}(:,ii)-B0{mdl,2}(:,ii))
    end
end

figure(23)
t = tiledlayout(3,2);
for ii = 1:3
    for mdl = 1:length(maxsz)
        nexttile(tilenum(t,ii,1))
        hold on
        plot(Be{mdl,1}(:,ii)-B0{mdl,1}(:,ii))
        nexttile(tilenum(t,ii,2))
        hold on
        plot(Be{mdl,2}(:,ii)-B0{mdl,2}(:,ii))
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
reference = B0{end,2};
norm_diff_B = zeros(size(V0,1), size(V0,2), 3);
for mdl = 1:size(B0,1)
    for ll = 1:size(V0,2)
        norm_diff_B(mdl,ll,:) = vecnorm(B0{mdl,ll}-reference,2,1)./vecnorm(reference,2,1);
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

for mdl = 1:size(B0,1)
    for ll = 1:size(B0,2)
        norm_diff_BB0(mdl,ll,:) = vecnorm(B0{mdl,ll}-B{mdl,ll},2,1)./vecnorm(B0{mdl,ll},2,1);
        norm_diff_B0Be(mdl,ll,:) = vecnorm(B0{mdl,ll}-Be{mdl,ll},2,1)./vecnorm(B0{mdl,ll},2,1);
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



