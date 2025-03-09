classdef DiscreteWaveletBasis < matlab.mixin.Copyable
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        coeff_matrix
        % basis
        % basis_vecs
        % elem_centers
        max_orders
        wname
        Avecs
        Dvecs
        n_elems
    end

    methods
        function obj = DiscreteWaveletBasis(max_orders)
            obj.wname = 'haar';
            % n = [8,8,8];
            obj.max_orders = max_orders;
            
            % provisory
            obj.coeff_matrix = obj.calc_DWT_coeffs(max_orders);
        end

        function calc_approx_detail_basis(obj)
            obj.Avecs = cell(length(obj.max_orders),1);
            obj.Dvecs = cell(length(obj.max_orders),1);
            for in = 1:length(obj.max_orders)
                nn = obj.max_orders(in);
                obj.Avecs{in} = obj.calc_approx_vec(nn);
                obj.Dvecs{in}  = obj.calc_detailed_vec(nn);
            end
        end

        function [basis, coeffs_matrix] =  make_basis(obj, elem_centers)
            obj.n_elems = size(elem_centers,1);
            % names = {'r', 'theta', 'z'};
            on_vec{1} = linspace(0,1,2*obj.max_orders(1));
            on_vec{2} = linspace(-pi,pi,2*obj.max_orders(2));
            on_vec{3} = linspace(-1,1,2*obj.max_orders(3));
            basis_vecs = cell(length(obj.max_orders),1);

            for in = 1:length(obj.max_orders)
                basis_vecs{in} = obj.interp_on_elem_centers_vec( ...
                    elem_centers(:,in), ...
                    on_vec{in}, obj.Avecs{in}, obj.Dvecs{in});
            end
            
            basis = obj.combine_basis(basis_vecs);
            [basis, nonzeros_idx] = obj.remove_zeros_basis(basis);
            % update coeffs_matrix
            coeffs_matrix = obj.coeff_matrix(nonzeros_idx,:);
            obj.coeff_matrix = coeffs_matrix;
        end



        function basis_vec = interp_on_elem_centers_vec(obj, elem_centers_vec, on_vec, Avec, Dvec)


            %%
            % [new_elem_centers] = helpers.cylindrical_elem_centers(elem_centers);
            % 
            % x = linspace(0,1,2*n(1));
            % y = linspace(-pi,pi,2*n(2));
            % z = linspace(-1,1,2*n(3));
            
            method = 'linear';
            interpolation_fun = @(x,v,xq) interp1(x, v, xq, method);
            
            n = size(Avec,1);
            % x
            vA = zeros(obj.n_elems,n);
            vD = zeros(obj.n_elems,n);
            for ii = 1:n
                vA(:,ii) = interpolation_fun(on_vec,Avec(ii,:),elem_centers_vec);
                vD(:,ii) = interpolation_fun(on_vec,Dvec(ii,:),elem_centers_vec);
            end
            
            basis_vec = [vA, vD];
        end

        function [coeffs_matrix] = calc_DWT_coeffs(obj, n)
            coeffs_matrix = zeros(prod(n),3);
            
            idx = 0;
            for ix = 1:2*n(1)
                for iy = 1:2*n(2)
                    for iz = 1:2*n(3)
                        idx = idx+1;
                        coeffs_matrix(idx,:) = [ix, iy, iz];
                    end
                end
            end
        end

        function Avec = calc_approx_vec(obj, nn)
            A = zeros(1,nn);
            D = zeros(1,nn);
            Avec = zeros(nn,2*nn);
        
            for ii = 1:nn
                a = A;
                d = D;
                a(ii) = 1;
                Avec(ii,:) = idwt(a,d,obj.wname);
        
            end
        end

        function Dvec = calc_detailed_vec(obj, nn)    
            A = zeros(1,nn);
            D = zeros(1,nn);
            Dvec = zeros(nn,2*nn);
        
            for ii = 1:nn
                a = A;
                d = D;
                d(ii) = 1;
                Dvec(ii,:) = idwt(a,d,obj.wname);
            end
        end

        function basis = combine_basis(obj, basis_vecs)
            basis_x = basis_vecs{1};
            basis_y = basis_vecs{2};
            basis_z = basis_vecs{3};
            
            basis = zeros(obj.n_elems,prod(2*obj.max_orders));
            
            idx = 0;
            for ix = 1:2*obj.max_orders(1)
                for iy = 1:2*obj.max_orders(2)
                    for iz = 1:2*obj.max_orders(3)
                        idx = idx+1;
                        basis(:,idx) = basis_x(:,ix).*basis_y(:,iy).*basis_z(:,iz);
                    end
                end
            end
        end

        function [basis, nonzeros_idx] = remove_zeros_basis(obj, basis)
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
end