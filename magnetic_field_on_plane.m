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
phantom.max_el_sz = 0.005;
phantom.maxsz = 0.01;
phantom.background = 0.503;

current_ampl = 10e-3;
eit = EIT(phantom, current_ampl);
eit.calc_elem_current();

integral_values = [];

%%
current_model.fmdl = eit.fwd_model;
current_model.elem_centers = eit.elem_centers;
current_model.elem_volumes = eit.elem_volumes;
current_model.elem_curr = eit.elem_curr(:,:,1);

coil_info.center = [phantom.radius + 0.01, 0, phantom.elec_vert_position];
coil_info.radius = 0.005;
coil_info.orientation = [0,pi/2, 0];

integral_values(:,1) = calc_coils_flux(coil_info, current_model);

%%
coil_info.radius = 0.005;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [pi/2,0, 0];

integral_values(:,2) = calc_coils_flux(coil_info, current_model);

%%
coil_info.radius = 0.005;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [0, 0, pi/2];

integral_values(:,3) = calc_coils_flux(coil_info, current_model);

%%

plot(integral_values)

%%


function integral_values = calc_coils_flux(coil_info, current_model)
    center = coil_info.center;
    radius = coil_info.radius;
    orientation = coil_info.orientation;
    
    mother_coil = Coil(center, radius, orientation);
    coils = CoilSystem(mother_coil, 16);
    
    
    % show model with coils
    coils.show();
    
    hold on
    show_fem(current_model.fmdl)
    hold off
    
    %%
    % current_model.elem_centers = elem_centers;
    % current_model.e_curr = e_curr;
    % current_model.elem_volumes = elem_volumes;
    
    integral_values = coils.calc_coil_integrals(current_model,1);
end

% 
%     strct = struct('B_positions_struct', [], 'B', [], 'abs_B', []);
%     B_results.vert_plane_tangential = strct;
%     B_results.horz_plane_tangential = strct;
%     B_results.vert_plane_sagital = strct;
% 
%     n_points = 25;
%     if ii ~= 1
%         what_2_plot.fem = false;
%     %     what_2_plot.abs_B = true;
%     %     what_2_plot.B = true;
%     end
% 
%     %% vert plane tangential to phantom
%     B_positions_struct = calc_position_vertical_plane_tangential_to_phantom(n_points, phantom);
%     B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
%     abs_B = vecnorm(B,2,2);
% 
%     plots_done{ii,1} = make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot);
% 
%     B_results.vert_plane_tangential.B_positions_struct = B_positions_struct;
%     B_results.vert_plane_tangential.B = B;
%     B_results.vert_plane_tangential.abs_B = abs_B;
% 
%     %% horz plane tangential to phantom
%     B_positions_struct = calc_position_horz_plane_tangential_to_phantom(n_points, phantom);
%     B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
%     abs_B = vecnorm(B,2,2);
% 
%     plots_done{ii,2} = make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot);
% 
%     B_results.horz_plane_tangential.B_positions_struct = B_positions_struct;
%     B_results.horz_plane_tangential.B = B;
%     B_results.horz_plane_tangential.abs_B = abs_B;
% 
%     %% vert plane sagital to phantom
%     [B_positions_struct] = calc_position_vert_plane_sagital_to_phantom(n_points, phantom);
%     B = helpers.calc_B_at_points(B_positions_struct.B_positions, elem_centers, e_curr, elem_volumes);
%     abs_B = vecnorm(B,2,2);
% 
%     plots_done{ii,3} = make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot);
% 
%     B_results.vert_plane_sagital.B_positions_struct = B_positions_struct;
%     B_results.vert_plane_sagital.B = B;
%     B_results.vert_plane_sagital.abs_B = abs_B;
% 
%     %%
%     all_results{ii,1} = max_el_sz;
%     all_results{ii,2} = B_results;
% 
% end
% %%
% save('FEM_B_convergency', 'all_results', 'phantom', 'max_el_sz_list')
% 
% plane_str = {'vert_plane_tangential', 'horz_plane_tangential', 'vert_plane_sagital'};
% 
% fun_strct = @(c,plane,what) c.(plane).(what);
% what_2_plot.fem = false;
% what_2_plot.abs_B = false;
% what_2_plot.B = true;
% 
% for ii = 1:length(max_el_sz_list)
%     B_results = all_results{ii,2};
%     for str_ = plane_str
%         B_positions_struct = fun_strct(B_results, str_{1}, 'B_positions_struct');
%         B = fun_strct(B_results, str_{1}, 'B');
%         make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot);
%     end
% 
% end
% 
% %%
% 
% plane_normals =[0,1,0; ...
%                 0,0,1; ...
%                 1,0,0];
% 
% 
% %%
% fun_strct = @(c, str) [c.vert_plane_tangential.(str), c.horz_plane_tangential.(str), c.vert_plane_sagital.(str)];
% 
% values_diff = zeros(length(max_el_sz_list), 3);
% norm_values = zeros(length(max_el_sz_list), 3);
% values = zeros(length(max_el_sz_list), 3);
% 
% % last maxsz is the reference
% for ii = 1:length(max_el_sz_list)
%     values_diff(ii,:) = vecnorm(fun_strct(all_results{ii,2}, 'abs_B')-fun_strct(all_results{end,2}, 'abs_B'),2,1);
%     norm_values(ii,:) = vecnorm(fun_strct(all_results{ii,2}, 'abs_B')-fun_strct(all_results{end,2}, 'abs_B'),2,1)./vecnorm(fun_strct(all_results{end,2}, 'abs_B'),2,1);
%     values(ii,:) = vecnorm(fun_strct(all_results{ii,2}, 'abs_B'),2,1);
% end
% 
% % for ii = 1:length(max_el_sz_list)
% %     all_B = fun_strct(all_results{ii,2}, 'B');
% %     fluxes(ii,:) = sum([all_B(:,1:3)*plane_normals(1,:)', ...
% %                         all_B(:,4:6)*plane_normals(2,:)', ...
% %                         all_B(:,7:9)*plane_normals(3,:)'],1);
% % end
% 
% figure(50)
% clf
% hold on
% plot(max_el_sz_list, values)
% plot(max_el_sz_list, values, '.')
% set(gca, 'xscale', 'log', 'yscale', 'log')
% xlabel('max electrode size / m')
% ylabel('Norm. values / T')
% 
% %%
% 
% x = log10(max_el_sz_list);
% y = log10(mean(values_diff,2));
% xx = log10(logspace(-5,max(x),500));
% p = polyfit(x,y,1);
% yy = polyval(p,xx);
% 
% figure(100)
% clf
% hold on
% plot(max_el_sz_list, values_diff, 'DisplayName', 'B_1')
% plot(max_el_sz_list, values_diff, 's')
% plot(10.^xx,10.^yy,'r','DisplayName','fit')
% set(gca, 'xscale', 'log', 'yscale', 'log')
% 
% ylabel('||B-B_0|| / ||B_0||')
% xlabel('h^{el}_{max} / m')
% legend()
% 
% x_ = 10.^mean(log10(xlim));
% y_ = 10.^mean(log10(ylim));
% text(x_,y_,sprintf('at h=%.2g m ||B-B_0||/||B_0||==0.01', 10^((log10(0.01)-p(2))/p(1))))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%% function declarations


% function B_position_plot3(B_positions)
% plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), 'ro','MarkerFaceColor','r');
% end
% 
% function [B_positions_struct] = calc_position_vertical_plane_tangential_to_phantom(n_points, phantom)
% x_dir = linspace(-phantom.phantom_radius, phantom.phantom_radius, n_points);
% y_dir = phantom.phantom_radius;
% z_dir = linspace(0, phantom.phantom_height, n_points);
% [xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
% B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];
% 
% B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);
% 
% [xx, zz] = meshgrid(x_dir, z_dir);
% 
% B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
% B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
% end
% 
% function [B_positions_struct] = calc_position_horz_plane_tangential_to_phantom(n_points, phantom)
% x_dir = linspace(-phantom.phantom_radius, phantom.phantom_radius, n_points);
% y_dir = linspace(phantom.phantom_radius, 3*phantom.phantom_radius, n_points);
% z_dir = phantom.phantom_height/2;
% 
% [xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
% B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];
% 
% B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);
% 
% [xx, zz] = meshgrid(x_dir, y_dir);
% 
% B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
% B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
% end
% 
% function [B_positions_struct] = calc_position_vert_plane_sagital_to_phantom(n_points, phantom)
% x_dir = 0;
% y_dir = linspace(phantom.phantom_radius, 3*phantom.phantom_radius, n_points);
% z_dir = linspace(0, phantom.phantom_height, n_points);
% 
% [xx, yy, zz] = meshgrid(x_dir, y_dir, z_dir);
% B_positions_struct.B_positions = [xx(:), yy(:), zz(:)];
% 
% B_positions_struct.plot3 = @() B_position_plot3(B_positions_struct.B_positions);
% 
% [xx, zz] = meshgrid(z_dir, y_dir);
% 
% B_positions_struct.surf = @(abs_B) surf(xx,zz,reshape(abs_B,size(xx)));
% B_positions_struct.contourf = @(Bi) contourf(xx, zz, reshape(Bi,size(xx)));
% end
% 
% function  [plots_done] = make_plots(fmdl, B_positions_struct, abs_B, B, what_2_plot)
% plots_done = {[],[],[]};
% if what_2_plot.fem
%     figure()
%     clf
%     hold on
%     show_fem(fmdl,[0,1.012])
%     B_positions_struct.plot3()
%     hold off
%     plots_done{1} = gcf;
% end
% 
% if what_2_plot.abs_B
%     figure()
%     % B_positions_struct.surf(abs_B)
%     B_positions_struct.contourf(abs_B); colorbar();
%     plots_done{2} = gcf;
% end
% 
% if what_2_plot.B
%     figure()
%     tiledlayout(1,3)
%     for ii = 1:3
%         nexttile
%         B_positions_struct.contourf(B(:,ii))
%         colorbar()
%     end
%     plots_done{3} = gcf;
% end
% end