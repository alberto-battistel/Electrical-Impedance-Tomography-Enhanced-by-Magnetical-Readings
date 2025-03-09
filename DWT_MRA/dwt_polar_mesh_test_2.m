
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

wname = 'haar';
n = 8;
A = zeros(1,n);
D = zeros(1,n);
Avec = zeros(n,2*n);
Dvec = zeros(n,2*n);

figure(1)
tiledlayout
clf
nexttile
hold on
for ii = 1:n
    a = A;
    d = D;
    a(ii) = 1;
    Avec(ii,:) = idwt(a,d,wname);
    plot(Avec(ii,:))
end

nexttile
hold on
for ii = 1:n
    a = A;
    d = D;
    d(ii) = 1;
    Dvec(ii,:) = idwt(a,d,wname);
    plot(Dvec(ii,:))
end

[new_elem_centers] = helpers.cylindrical_elem_centers(elem_centers, '');

x = linspace(-1,1,2*n);
y = linspace(0,pi,2*n);

method = 'linear'
% interpolation_fun = @interp1;
interpolation_fun = @(x,v,xq) interp1(x, v, xq, method);
% x
vAx = zeros(length(elem_centers),n);
for ii = 1:n
    vAx(:,ii) = interpolation_fun(x,Avec(ii,:),new_elem_centers(:,1));
end

% figure(2)
% clf
% tiledlayout()
% for ii = 1:n
%     img.elem_data = vAx(:,ii);
%     nexttile
%     show_fem(img, [1])
% end

vDx = zeros(length(elem_centers),n);
for ii = 1:n
    vDx(:,ii) = interpolation_fun(x,Dvec(ii,:),new_elem_centers(:,1));
end

% figure(3)
% clf
% tiledlayout()
% for ii = 1:n
%     img.elem_data = vDx(:,ii);
%     nexttile
%     show_fem(img, [1])
% end

% y
vAy = zeros(length(elem_centers),n);
for ii = 1:n
    vAy(:,ii) = interpolation_fun(y,Avec(ii,:),new_elem_centers(:,2));
end

% figure(4)
% clf
% tiledlayout()
% for ii = 1:n
%     img.elem_data = vAy(:,ii);
%     nexttile
%     show_fem(img, [1])
% end

vDy = zeros(length(elem_centers),n);
for ii = 1:n
    vDy(:,ii) = interpolation_fun(y,Dvec(ii,:),new_elem_centers(:,2));
end

% figure(5)
% clf
% tiledlayout()
% for ii = 1:n
%     img.elem_data = vDy(:,ii);
%     nexttile
%     show_fem(img, [1])
% end


%% tensor product-like basis

basis_x = [vAx,vDx];
basis_y = [vAy,vDy];

basis = zeros(length(elem_centers),n*n);
coeffs_matrix = zeros(n*n,2);

idx = 0;
for ix = 1:2*n
    for iy = 1:2*n
        idx = idx+1;
        basis(:,idx) = basis_x(:,ix).*basis_y(:,iy);
        coeffs_matrix(idx,:) = [ix, iy];
    end
end

% some basis may be zeros
l = zeros(1,size(basis,2));
for ii = 1:size(basis,2)
    l(ii) = ~all(basis(:,ii) == 0);
end
l = logical(l);

basis = basis(:,l);

fprintf('Basis size = %d, ~(2n * 2n)\n', size(basis,2))



%% plot some basis
figure(5)
clf
tiledlayout()
id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
warning('off',id)

for ii = [1:32]
    img.elem_data = basis(:,ii);
    nexttile
    show_fem(img, [1])
end
warning('on',id)

figure(6)
clf
tiledlayout()
id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
warning('off',id)

for ii = [1:16:256]
    img.elem_data = basis(:,ii);
    nexttile
    show_fem(img, [1])
end
warning('on',id)

%%
sol = zeros(length(v),2);
names = {};

%% least square, 
% not good
names{1} = 'Least square';
x1 = basis\v;
sol(:,1) = basis*x1;
fprintf('least squares\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,1)), n)

%% inner product, numerically stable, 
% but not good
x2 = sum(v.*basis,1);
sol(:,2) = basis*x2';
fprintf('inner product\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,2)), n)

%% lasso,
% good
names{3} = 'Lasso';

[x3,FitInfo] = lasso(basis,v);

sol(:,3) = basis*x3(:,10);
fprintf('lasso\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,3)), n)

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
