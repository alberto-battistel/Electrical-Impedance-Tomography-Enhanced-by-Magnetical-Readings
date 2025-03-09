
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
% wname = 'db4'
n = 10;
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

[new_elem_centers] = helpers.shift_elem_centers(elem_centers, [1,1], [2*n,2*n]);
x = 1:2*n;
y = 1:2*n;

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


%%
sol = zeros(length(v),2);

%%

A = [vAx,vDx,vAy,vDy];

x1 = A\v;
x2 = sum(v.*A,1);

sol(:,1) = A*x1;
sol(:,2) = A*x2';
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,1)), n)
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,2)), n)

%%

basis_x = [vAx,vDx];
basis_y = [vAy,vDy];

basis = zeros(length(elem_centers),n*n);

idx = 0;
for ix = 1:2*n
    for iy = 1:2*n
        idx = idx+1;
        basis(:,idx) = basis_x(:,ix).*basis_y(:,iy);
    end
end

%%
x3 = basis\v;
sol(:,3) = basis*x3;
fprintf('least squares\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,3)), n)

%% inner product

x4 = sum(v.*basis,1);
sol(:,4) = basis*x4';
fprintf('inner product\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,4)), n)

%%
[U,S,V] = svd(basis);
l = find(diag(S)>=1e-8*S(1,1),1,'last')
bb = U(:,1:2*n.^2);
% bb = U(:,1:l);

x5 = sum(v.*bb,1);
sol(:,5) = bb*x5';
fprintf('inner product svd\n')
fprintf('L2= %f with n=%d\n', vecnorm(v-sol(:,5)), n)

x6 = bb\v;
sol(:,6) = bb*x6;
fprintf('least squares svd\n')
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