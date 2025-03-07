home
clear
close all
init_eidors()

addpath('..')
%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 0.1;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001; %0.0025
phantom.maxsz = 0.0025;
phantom.background = 0.503; % muscle at 1 MHz
phantom.extra = {'', ''};
phantom.extra_format = {'', ''};

n_coils = 1;
conductivity = phantom.background;

current_ampl = 10e-3;
freq = 1e6;

max_el_sz = 0.001:0.001:0.003;


models = struct('eit', [], 'max_el_sz', [], 'coil_detectors', []);

tic
for ii = 1:length(max_el_sz)
    phantom.max_el_sz = max_el_sz(ii);

    fprintf('phantom max_el_sz: %.4f m\n', phantom.max_el_sz)
    
    models(ii).max_el_sz = phantom.max_el_sz;

    models(ii).eit = EIT(phantom, current_ampl);
end


%%

for model = 1:length(models)
    models(model).eit.calc_elem_current(model);
end


%%
B_position = [0-0.005, ...
    models(1).eit.phantom.radius + 0.01-0.005, ...
    models(1).eit.phantom.elec_vert_position-0.005];


figure(1)
clf
models(1).eit.show_fem
hold on
plot3(B_position(1), B_position(2), B_position(3), 'sm')
xlabel('x')
ylabel('y')
zlabel('z')

%%
B_fem = zeros(length(models),3);
for model = 1:length(models)
    elem_centers = models(model).eit.elem_centers;
    elem_currents = models(model).eit.elem_currents;
    elem_volumes = models(model).eit.elem_volumes;
    B_fem(model,:) = helpers.calc_B_at_points(B_position, elem_centers, elem_currents, elem_volumes);

end

%%
figure(543)
clf
tiledlayout(3,1)
for comp = 1:3
    nexttile
    plot(max_el_sz, B_fem(:,comp))
end




%%
N2 = 251;
[xx,yy,zz] = meshgrid(linspace(-phantom.radius,phantom.radius,N2), ...
                    linspace(-phantom.radius,phantom.radius,N2), ...
                    linspace(0,phantom.height,N2));


%%

[interp_voltages] = zeros(length(xx),length(yy),length(zz),length(models));

for model = 1:length(models)
    elem_nodes = models(model).eit.fwd_model.nodes;
    elem_values = models(model).eit.volt_strct.volt(:,1);

    interp_voltages(:,:,:,model) = interpolate_elem_values(xx, yy, zz, elem_nodes, elem_values);
end

rel_diff_interp_voltages = log10(abs(interp_voltages./interp_voltages(:,:,:,1)-1));

%%
figure(97634)
clf
tiledlayout(1,length(models))

z_idx = floor(N2/2);

min_val = min(interp_voltages(:,:,z_idx,:),[],'all');
max_val = max(interp_voltages(:,:,z_idx,:),[],'all');
levels = linspace(min_val, max_val, 11);
for model = 1:length(models)
    nexttile

    contourf(squeeze(interp_voltages(:,:,50,model)))
    colorbar
    clim([min_val max_val])
end

%%
figure(2598)
clf
tiledlayout(1,length(models)-1)

z_idx = 50;
for model = 2:length(models)
    nexttile
    histogram(rel_diff_interp_voltages(:,:,:,model),51)
    xlim([-6,0])
end

%% electric field
% dx = diff(xx(1,1:2,1));
% dy = diff(yy(1:2,1,1));
% dz = diff(zz(1,1,1:2));
% [ex,ey,ez] = gradient(interp_voltages(:,:,:,1), dx,dy,dz);
% 
% ix = conductivity*ex;
% iy = conductivity*ey;
% iz = conductivity*ez;

%%


elem_volumes = dx*dy*dz;

% elem_currents = [ix(:),iy(:),iz(:)];

elem_centers = [xx(:),yy(:),zz(:)];

B_fd = zeros(length(models),3);
for model = 1:length(models)

    [elem_currents] = calc_current_from_voltages(interp_voltages(:,:,:,model), dx,dy,dz, conductivity);

    B_fd(model,:) = helpers.calc_B_at_points(B_position, elem_centers, elem_currents, elem_volumes);
end

%%
figure(544+1)
clf
tiledlayout(3,1)
for comp = 1:3
    nexttile
    plot(max_el_sz, B_fd(:,comp))
end

%%
[interp_volumes] = zeros(length(xx),length(xx),length(models));

for model = 1:length(models)
    elem_centers = models(model).eit.elem_centers;
    elem_values = models(model).eit.elem_volumes;

    interp_volumes(:,:,model) = interpolate_elem_values(xx, yy, zz, elem_centers, elem_values);

end

rel_diff_interp_volumes = log10(abs(interp_volumes./interp_volumes(:,:,1)-1));

%%
[interp_dB] = zeros(length(xx),length(xx),3,length(models));
point = B_position;


for model = 1:length(models)
    elem_centers = models(model).eit.elem_centers;
    r = point - elem_centers;
    j_times_vol = models(model).eit.elem_currents.*models(model).eit.elem_volumes;
    dB = cross(j_times_vol, r)./vecnorm(r,2,2).^3; % magnetic flux density, in T
    interp_dB(:,:,:,model) = interpolate_elem_values(xx, yy, zz, elem_centers, dB);

end

rel_diff_interp_dB = log10(abs(interp_dB./interp_dB(:,:,:,1)-1));

%%
figure(563453)
clf
tiledlayout(3,length(models)-1)

for comp = 1:3
    for model = 2:length(models)
        nexttile
        histogram(rel_diff_interp_dB(:,:,comp,model),21)
        xlim([-4,4])
    end
end

%%
figure(4534)
clf
tiledlayout(3,length(models))

for comp = 1:3
    min_val = min(interp_dB(:,:,comp,:),[],'all');
    max_val = max(interp_dB(:,:,comp,:),[],'all');
    levels = linspace(min_val, max_val, 11);
    for model = 1:length(models)
        nexttile
        contourf(interp_dB(:,:,comp,model),levels)
        clim([min_val max_val])
        xlim(N2/2+[-50,+50])
        ylim(N2.*[0.9,1])
        colorbar
    end
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


