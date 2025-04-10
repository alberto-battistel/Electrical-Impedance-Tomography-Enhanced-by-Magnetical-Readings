home
clear
close all
init_eidors()

addpath('..')
addpath('../progressbar/')

%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 1.*phantom.radius;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.005;
phantom.maxsz = 0.01; % final one to use 0.001
phantom.background = 0.503; % muscle at 1 MHz

coil_radia = [0.005, 0.01, 0.015, 0.02, 0.025];


current_ampl = 10e-3;
freq = 1e6;
eit = EIT(phantom, current_ampl);

eit.show_fem()
eit.calc_elem_current();



%%
current_model.fmdl = eit.fwd_model;
current_model.img = eit.img;
current_model.elem_centers = eit.elem_centers;
current_model.elem_volumes = eit.elem_volumes;
current_model.elem_currents = eit.elem_currents;
current_model.measurement_idx = 1:phantom.n_elec;

integral_values = zeros(phantom.n_elec, phantom.n_elec, 3, length(coil_radia));

%%
progressbar
for ii = 1:length(coil_radia)
    %% horizontal tangential coils
    integral_values(:,:,1,ii) = integral_on_horizontal_tangential_coils(phantom, coil_radia(ii), current_model);
    
    %% vertical perpendicolar coils
    integral_values(:,:,2,ii) = integral_on_vertical_perpendicolar_coils(phantom, coil_radia(ii), current_model);
    
    %% horizontal coplanar coils
    integral_values(:,:,3,ii) = integral_on_horizontal_coplanar_coils(phantom, coil_radia(ii), current_model);

    progressbar(ii/length(coil_radia))
end

%%

magnetic_voltages = 2*pi*freq*reshape(integral_values, [], 3, length(coil_radia));


%% plot rms voltage coils vs coil radius

titles = {'Tangential Coils', 'Perpendicolar Coils', 'Coplanar Coils'};

figure(1256)
tiledlayout(3,1)
for ii = 1:3
    nexttile
    bar(coil_radia*100, rms(squeeze(magnetic_voltages(:,ii,:))))
    % set(gca, 'yscale', 'log')
    xlim([0, 3])
    xlabel('Coil Radium / cm')
    ylabel('Coil Voltage / V_{rms}')
    title(titles{ii})
end

%% plot the coil voltage for the first coil

figure(167456)
tiledlayout(3,1)
for ii = 1:3
    nexttile
    plot(1:1:phantom.n_elec^2, magnetic_voltages(:,ii,1))
    xlim([0, 257])
    xlabel('Coil Radium / cm')
    ylabel('Coil Voltage / V')
    title(titles{ii})
end



%% function declarations


function integral_values = calc_coils_flux(coil_info, current_model)
    center = coil_info.center;
    radius = coil_info.radius;
    orientation = coil_info.orientation;
    
    mother_coil = Coil(center, radius, orientation);
    coils = CoilSystem(mother_coil, length(current_model.measurement_idx));
    
    
    % show model with coils
    % coils.show();
    
    % hold on
    % show_fem(current_model.img)
    % hold off
        
    integral_values = coils.calc_coil_integrals(current_model,current_model.measurement_idx);
end


function integral_values = integral_on_horizontal_tangential_coils(phantom, coil_radius, current_model)
%% horizontal tangential coils
coil_info.center = [phantom.radius + 0.01, 0, phantom.elec_vert_position];
coil_info.radius = coil_radius;
coil_info.orientation = [0,pi/2, 0];

integral_values = calc_coils_flux(coil_info, current_model);
end


function integral_values = integral_on_vertical_perpendicolar_coils(phantom, coil_radius, current_model)
%% vertical perpendicolar coils
coil_info.radius = coil_radius;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [pi/2,0, 0];

integral_values = calc_coils_flux(coil_info, current_model);
end


function integral_values = integral_on_horizontal_coplanar_coils(phantom, coil_radius, current_model)
%% horizontal coplanar coils
coil_info.radius = coil_radius;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [0, 0, pi/2];

integral_values = calc_coils_flux(coil_info, current_model);
end

