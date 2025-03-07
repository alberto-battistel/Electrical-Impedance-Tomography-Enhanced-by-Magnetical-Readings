
close all
addpath('..')
mk_model_model_str = 'd2c';

imdl = mk_common_model(mk_model_model_str,16);

img = mk_image(imdl);

elem_centers = interp_mesh(img.fwd_model, 0); % center of elements

%%

pp = [.1,0.2];
radius = 0.3;

l = sum((elem_centers-pp).^2,2) <= radius.^2;

v = zeros(length(elem_centers),1);
v(l) = 2.5;

%%
n = 16;
fun = 1.1*ones(1,n);
fun(7:13) = -2.5;

[c,l] = wavedec(fun,3,"haar");

basis_ = zeros(length(c), size(fun,2));

C = zeros(size(c));
for ii = 1:length(c)
    cc = C;
    cc(ii) = 1;
    basis_(ii,:) = waverec(cc,l,"haar");
end

figure(345)
clf
tiledlayout()

l1 = 1;
for ii = l(1:end-1)
    nexttile
    l2 = l1+ii-1;
    plot(basis_(l1:l2,:)')
    l1 = l2+1;
    title(num2str(ii))
end

%% Interpolate wavelets on mesh
[new_elem_centers] = helpers.cylindrical_elem_centers(elem_centers);

x = linspace(0,1,n);
y = linspace(-pi,pi,n);

method = 'linear'
interpolation_fun = @(x,v,xq) interp1(x, v, xq, method);

% x
basis_x = zeros(length(elem_centers),n);
for ii = 1:n
    basis_x(:,ii) = interpolation_fun(x,basis_(ii,:),new_elem_centers(:,1));
end

% y
basis_y = zeros(length(elem_centers),n);
for ii = 1:n
    basis_y(:,ii) = interpolation_fun(y,basis_(ii,:),new_elem_centers(:,2));
end

% tensor product-like basis
basis = zeros(length(elem_centers),n*n);

idx = 0;
for ix = 1:n
    for iy = 1:n
        idx = idx+1;
        basis(:,idx) = basis_x(:,ix).*basis_y(:,iy);
    end
end

fprintf('Basis size = %d, (n*n)\n', size(basis,2))

%% plot some basis
figure(5)
clf
tiledlayout()
id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
warning('off',id)

for ii = [1,5,25,50,72,125,181,233]
    img.elem_data = basis(:,ii);
    nexttile
    show_fem(img, [1])
end
warning('on',id)

%%
sol = zeros(length(v),2);
names = {};

%% least square, 
% it gives numerical problems
names{1} = 'Least square';
x1 = basis\v;
sol(:,1) = basis*x1;
fprintf('least squares\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,1)), n)

%% inner product, numerically stable, 
% but not good
names{2} = 'Inner product';
x2 = sum(v.*basis,1);
sol(:,2) = basis*x2';
fprintf('inner product\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,2)), n)

%% proper orthogonal reduction?, 
% good
names{3} = 'POR?';
[coeff,score,latent,tsquared,explained,mu] = pca(basis);

l = find(cumsum(explained)>=99.5,1,'first')
ss = zeros(size(score));
ss(:,1:l) = score(:,1:l);
% A0 = score*coeff'+mu;
aa = score+mu;
A = ss*coeff'+mu;

x3 = A\v;
sol(:,3) = A*x3;
fprintf('proper orthogonal reduction?\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,3)), n)

%% proper orthogonal reduction? with new basis, 
% good but with waves
names{4} = 'POR? new basis';
l = find(cumsum(explained)>=99.9,1,'first')

A = score(:,1:l)+mu(1:l);

x4 = A\v;
sol(:,4) = A*x4;
fprintf('proper orthogonal reduction? with new basis\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,4)), n)

%% lasso,
% good
names{5} = 'Lasso';

[x5,FitInfo] = lasso(basis,v);

sol(:,5) = basis*x5(:,10);
fprintf('lasso\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,5)), n)

%% lasso with proper orthogonal reduction? new basis, 
% good
names{6} = 'Lasso with POR new basis';
[x6,FitInfo] = lasso(A,v);

sol(:,6) = A*x6(:,10);
fprintf('lasso properly reduced\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,6)), n)


%%


figure(46)
clf
tiledlayout(size(sol,2),4)

id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
warning('off',id)

for ii = 1:size(sol,2)
    img.elem_data = v;
    nexttile
    show_fem(img, [1])
    title(names{ii})

    img.elem_data = sol(:,ii);
    nexttile
    show_fem(img, [1])

    img.elem_data = v-sol(:,ii);
    nexttile
    show_fem(img, [1])

    img.elem_data = vecnorm(sol(:,ii)-v,2,2);
    nexttile
    show_fem(img, [1])

end

warning('on',id)