% from https://eidors3d.sourceforge.net/doc/index.html?eidors/graphics/matlab/show_current.html
% do_unit_test


% init_eidors
close all
clear

mu0 = 4*pi*1e-7;
n_electrodes = 8;

fmdl = mk_common_model('g2c2',n_electrodes);   % by defaul 10 A
fmdl.normalize_measurements = 0;
img = mk_image(fmdl); 
img.fwd_solve.get_all_meas = 1;

% 10 mA current
for ii = 1:length(img.fwd_model.stimulation)
    img.fwd_model.stimulation(ii).stim_pattern = 1e-3*img.fwd_model.stimulation(ii).stim_pattern;
end

elem_centers = interp_mesh(fmdl);
elem_areas = calc_area(fmdl.fwd_model.nodes, fmdl.fwd_model.elems);

angles = exp(1j*(pi/2-(2*pi/n_electrodes)*((0:n_electrodes-1)+0.5)));
coil_positions = 1.05*[real(angles)', imag(angles)'];

vh = fwd_solve(img);
B = zeros(n_electrodes,n_electrodes);
for i_volt = 1:n_electrodes
    e_curr = calc_elem_current(img, vh.volt(:,i_volt));

    for i_point = 1:length(coil_positions)
        point = coil_positions(i_point,:);
        r = point - elem_centers;
        % if need to check for zero division
        r_mag = norm(r);
    
        dB = mu0/(4*pi)*cross_prod_2d(e_curr.*elem_areas, r)./norm(r.^3); % magnetic flux density, in T
    
        B(i_point, i_volt) = sum(dB);
    end
end

figure(100)
tiledlayout(4,2)
for ii = 1:4
    nexttile
    show_current(img, vh.volt(:,ii));
    nexttile
    plot(B(:,ii))
end

%%


