function [basis, nonzeros_idx] = make_DWT_set(n, elem_centers)

wname = 'haar';
% n = [8,8,8];

n_elems = size(elem_centers,1);

% names = {'r', 'theta', 'z'};
on_vec{1} = linspace(0,1,2*n(1));
on_vec{2} = linspace(-pi,pi,2*n(2));
on_vec{3} = linspace(-1,1,2*n(3));


basis_vecs = cell(length(n),1);
for in = 1:length(n)
    nn = n(in);
    Avec = calc_approx_vec(nn);
    Dvec = calc_detailed_vec(nn);

    basis_vecs{in} = interp_on_elem_centers_vec(elem_centers(:,in), on_vec{in}, Avec, Dvec);
end

basis = combine_basis(basis_vecs);
[basis, nonzeros_idx] = remove_zeros_basis(basis);



function Avec = calc_approx_vec(nn)

% figure(1)
% clf
% tiledlayout(2,3)

% for in = 1:length(n)
    % nn = n(in);
    A = zeros(1,nn);
    D = zeros(1,nn);
    Avec = zeros(2*nn,2*nn);
    % nexttile
    % hold on
    for ii = 1:nn
        a = A;
        d = D;
        a(ii) = 1;
        Avec(ii,:) = idwt(a,d,wname);
        % plot(Avec(ii,:))
    end
    % title(names{in})
    % Avec{in} = Avec_;
% end
% legend('Approx.')
end

function Dvec = calc_detailed_vec(nn)
% for in = 1:length(n)
%     nn = n(in);
    A = zeros(1,nn);
    D = zeros(1,nn);
    Dvec = zeros(2*nn,2*nn);
    % nexttile
    % hold on
    for ii = 1:nn
        a = A;
        d = D;
        d(ii) = 1;
        Dvec(ii,:) = idwt(a,d,wname);
        % plot(Dvec(ii,:))
    end
    % Dvec{in} = Dvec;
% end
% legend('Detailed')
end

function basis_vec = interp_on_elem_centers_vec(elem_centers_vec, on_vec, Avec, Dvec)


%%
% [new_elem_centers] = helpers.cylindrical_elem_centers(elem_centers);
% 
% x = linspace(0,1,2*n(1));
% y = linspace(-pi,pi,2*n(2));
% z = linspace(-1,1,2*n(3));

method = 'linear';
interpolation_fun = @(x,v,xq) interp1(x, v, xq, method);

% x
vA = zeros(n_elems,length(Avec));
vD = zeros(n_elems,length(Dvec));
for ii = 1:n(1)
    vA(:,ii) = interpolation_fun(on_vec,Avec(ii,:),elem_centers_vec);
    vD(:,ii) = interpolation_fun(on_vec,Dvec(ii,:),elem_centers_vec);
end

basis_vec = [vA, vD];
end

%%
function basis = combine_basis(basis_vecs)
basis_x = basis_vecs{1};
basis_y = basis_vecs{2};
basis_z = basis_vecs{3};

basis = zeros(n_elems,prod(2*n));

idx = 0;
for ix = 1:2*n(1)
    for iy = 1:2*n(2)
        for iz = 1:2*n(3)
            idx = idx+1;
            basis(:,idx) = basis_x(:,ix).*basis_y(:,iy).*basis_z(:,iz);
        end
    end
end
end

function [basis, nonzeros_idx] = remove_zeros_basis(basis)
% some basis may be zeros
nonzeros_idx = zeros(1,size(basis,2));
for ii = 1:size(basis,2)
    nonzeros_idx(ii) = ~all(basis(:,ii) == 0);
end
nonzeros_idx = logical(nonzeros_idx);

basis = basis(:,nonzeros_idx);
% n_basis = size(basis,2);
% fprintf('Basis size = %d, ~(2n1 * 2n2 * 2n3)\n', n_basis)
end

end