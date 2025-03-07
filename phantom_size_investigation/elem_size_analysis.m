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
phantom.max_el_sz = []; %0.0025
phantom.maxsz = 0.005; % 0.0025 or 0.005
phantom.background = 0.503; % muscle at 1 MHz

% this 
phantom.extra_format = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
phantom.extra = phantom.extra_format(0.025, [0,0,phantom.elec_vert_position]);
phantom.ball = 0.503; % muscle at 1 MHz
% or this
% phantom.extra = {'', ''};
% phantom.extra_format = {'', ''};

n_coils = 1;

current_ampl = 10e-3;
freq = 1e6;

max_el_sz = 0.001:0.001:0.005;


models = struct('eit', [], 'max_el_sz', [], 'coil_detectors', []);

tic
for ii = length(max_el_sz):-1:1
    phantom.max_el_sz = max_el_sz(ii);

    fprintf('phantom max_el_sz: %.4f m\n', phantom.max_el_sz)
    
    models(ii).max_el_sz = phantom.max_el_sz;

    models(ii).eit = EIT(phantom, current_ampl);
end


%%

for model = 1:length(models)
    models(model).eit.calc_elem_current(1);
end


%%
% B_position = [0.05, models(1).eit.phantom.radius + 0.05, models(1).eit.phantom.elec_vert_position+0.05]; % B away from symmetry planes
B_position = [0, models(1).eit.phantom.radius + 0.01, models(1).eit.phantom.elec_vert_position]; % exactly on a symmetry plance

figure(1)
clf
models(1).eit.show_fem
hold on
plot3(B_position(1), B_position(2), B_position(3), 'sm')
xlabel('x')
ylabel('y')
zlabel('z')

%%
B = zeros(length(models),3);
for model = 1:length(models)
    elem_centers = models(model).eit.elem_centers;
    elem_currents = models(model).eit.elem_currents;
    elem_volumes = models(model).eit.elem_volumes;
    B(model,:) = helpers.calc_B_at_points(B_position, elem_centers, elem_currents, elem_volumes);

end
rel_diff_B = (B-B(1,:))./B(1,:)*100;

rel_norm_diff_B = vecnorm(B-B(1,:), 2,2)./vecnorm(B(1,:), 2,2);


%%

comps = {'x', 'y', 'z'};
figure(546)
clf
tiledlayout(2,3)
for comp = 1:3
    nexttile
    plot(max_el_sz, B(:,comp))
    title(comps{comp})
    ylabel('B / T')
end

for comp = 1:3
    nexttile
    plot(max_el_sz(1:end), rel_diff_B(1:end,comp))
    ylabel('relative error / %')
    xlabel('max_el_sz', 'Interpreter','none')
    title(comps{comp})
end

if ~isempty(phantom.extra{1})
    str = 'with ball';
else
    str = '';
end

sgtitle(sprintf('maxsz = %.5f %s', phantom.maxsz, str))


figure(5954)
clf
tiledlayout(2,3)
for comp = 1:3
    nexttile
    plot(max_el_sz, B(:,comp))
    title(comps{comp})
    ylabel('B / T')
end


    nexttile([1 3])
    plot(max_el_sz, rel_norm_diff_B)
    ylabel('rel. norm diff. B')
    xlabel('max_el_sz', 'Interpreter','none')
    


sgtitle(sprintf('maxsz = %.5f %s', phantom.maxsz, str))

%%

num_elems = zeros(length(models),1);
for ii = 1:length(models)
    num_elems(ii) = length(models(ii).eit.elem_centers);
end

figure(342516)

plot(max_el_sz, num_elems)
xlabel('max_el_sz', 'Interpreter','none')
ylabel('n. elements')
sgtitle(sprintf('maxsz = %.5f %s', phantom.maxsz, str))

%%
voltages = zeros(length(models),208);


for ii = 1:length(models)
    voltages(ii,:) = models(ii).eit.volt_strct.meas;
end

rel_diff_voltages = (voltages-voltages(1,:))./voltages(1,:)*100;

rel_norm_diff_voltages = vecnorm(voltages-voltages(1,:), 2,2)./vecnorm(voltages(1,:), 2,2);

figure(2516)

plot(max_el_sz, rel_norm_diff_voltages)
xlabel('max_el_sz', 'Interpreter','none')
ylabel('rel. norm diff. voltages')
sgtitle(sprintf('maxsz = %.5f %s', phantom.maxsz, str))


toc