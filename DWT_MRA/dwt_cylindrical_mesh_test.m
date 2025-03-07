
close all
addpath('..')


% cylinder r=0.1, h=0.2, c_base=[0,0,0] 
shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,0.1; 0.1); \n', ...
'solid bottom = plane(0,0,0;0,0,-1);\n' ...
'solid top    = plane(0,0,0.2;0,0,1);\n' ...
'solid mainobj= top and bottom and cyl -maxh=0.005;\n'];
elec_pos = [  0,  0,  0,   0,  0,  1];
elec_shape=[0,0.1,0.2];
elec_obj = {'bottom'};
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

figure(1)
show_fem(fmdl)

img = mk_image(fmdl,0);

elem_centers = interp_mesh(img.fwd_model, 0); % center of elements

%%

pp = [0.035,0.025,0.1];
radius = 0.025;

l = sum((elem_centers-pp).^2,2) <= radius.^2;

v = zeros(length(elem_centers),1);
v(l) = 2.5;

level = [inf,inf,0.1];
img.elem_data = v;
figure(8753)
show_slices(img, level)

%%

wname = 'haar';
n = [8,8,8];
names = {'r', 'theta', 'z'};

Avec = cell(3,1);
Dvec = cell(3,1);


figure(1)
clf
tiledlayout(2,3)

for in = 1:length(n)
    nn = n(in);
    A = zeros(1,nn);
    D = zeros(1,nn);
    Avec_ = zeros(2*nn,2*nn);
    nexttile
    hold on
    for ii = 1:nn
        a = A;
        d = D;
        a(ii) = 1;
        Avec_(ii,:) = idwt(a,d,wname);
        plot(Avec_(ii,:))
    end
    title(names{in})
    Avec{in} = Avec_;
end
legend('Approx.')

for in = 1:length(n)
    nn = n(in);
    A = zeros(1,nn);
    D = zeros(1,nn);
    Dvec_ = zeros(2*nn,2*nn);
    nexttile
    hold on
    for ii = 1:nn
        a = A;
        d = D;
        d(ii) = 1;
        Dvec_(ii,:) = idwt(a,d,wname);
        plot(Dvec_(ii,:))
    end
    Dvec{in} = Dvec_;
end
legend('Detailed')


% nexttile
% hold on
% for ii = 1:nn
%     a = A;
%     d = D;
%     d(ii) = 1;
%     Dvec(ii,:) = idwt(a,d,wname);
%     plot(Dvec(ii,:))
% end

%%
[new_elem_centers] = helpers.cylindrical_elem_centers(elem_centers);

x = linspace(0,1,2*n(1));
y = linspace(-pi,pi,2*n(2));
z = linspace(-1,1,2*n(3));

method = 'linear'

interpolation_fun = @(x,v,xq) interp1(x, v, xq, method);
% x
vAx = zeros(length(elem_centers),n(1));
vDx = zeros(length(elem_centers),n(1));
for ii = 1:n(1)
    vAx(:,ii) = interpolation_fun(x,Avec{1}(ii,:),new_elem_centers(:,1));
    vDx(:,ii) = interpolation_fun(x,Dvec{1}(ii,:),new_elem_centers(:,1));
end

% y
vAy = zeros(length(elem_centers),n(2));
vDy = zeros(length(elem_centers),n(2));
for ii = 1:n(2)
    vAy(:,ii) = interpolation_fun(y,Avec{2}(ii,:),new_elem_centers(:,2));
    vDy(:,ii) = interpolation_fun(y,Dvec{2}(ii,:),new_elem_centers(:,2));
end

% z
vAz = zeros(length(elem_centers),n(3));
vDz = zeros(length(elem_centers),n(3));
for ii = 1:n(3)
    vAz(:,ii) = interpolation_fun(z,Avec{3}(ii,:),new_elem_centers(:,3));
    vDz(:,ii) = interpolation_fun(z,Dvec{3}(ii,:),new_elem_centers(:,3));
end

%%

basis_x = [vAx,vDx];
basis_y = [vAy,vDy];
basis_z = [vAz,vDz];

basis = zeros(length(elem_centers),prod(2*n));

idx = 0;
for ix = 1:2*n(1)
    for iy = 1:2*n(2)
        for iz = 1:2*n(3)
            idx = idx+1;
            basis(:,idx) = basis_x(:,ix).*basis_y(:,iy).*basis_z(:,iz);
        end
    end
end

% some basis may be zeros
l = zeros(1,size(basis,2));
for ii = 1:size(basis,2)
    l(ii) = ~all(basis(:,ii) == 0);
end
l = logical(l);

basis = basis(:,l);
n_basis = size(basis,2)
fprintf('Basis size = %d, ~(2n1 * 2n2 * 2n3)\n', n_basis)


%% plot some basis
% to fix
levels = [inf,inf,0.1; ...
    inf, 0, inf; ...
    0 , inf, inf];

basis_to_see = 1:5;

figure(5)
clf
t1 = tiledlayout(3,length(basis_to_see),"TileSpacing", "none");
% id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
% warning('off',id)

for iii = 1:3
    for ii = 1:length(basis_to_see)
        img.elem_data = basis(:,ii);
        nexttile
    
        show_slices(img, levels(iii,:))
    end
end
% warning('on',id)

%%
sol = zeros(length(v),2);
names = {};

%% least square, 
% good
names{1} = 'Least square';
x1 = basis\v;
sol(:,1) = basis*x1;
fprintf('least squares\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,1)), n_basis)

%% lasso,
% good
names{2} = 'Lasso';

[x2,FitInfo] = lasso(basis,v);

sol(:,2) = basis*x2(:,1);
fprintf('lasso\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,2)), n_basis)

figure(623)
plot(FitInfo.MSE)

%%
figure(46)
clf
tiledlayout(size(sol,2),4)

level = [inf,inf,0.1];

for ii = 1:size(sol,2)
    img.elem_data = v;
    nexttile
    show_slices(img, level)
    title(names{ii})

    img.elem_data = sol(:,ii);
    nexttile
    show_slices(img, level)
    title('Reconstructed')

    img.elem_data = v-sol(:,ii);
    nexttile
    show_slices(img, level)
    title('Difference')

    img.elem_data = vecnorm(sol(:,ii)-v,2,2);
    nexttile
    show_slices(img, level)
    title('L2')

end


