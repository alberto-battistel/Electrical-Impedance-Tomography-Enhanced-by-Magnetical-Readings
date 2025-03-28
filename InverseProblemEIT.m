classdef (Abstract) InverseProblemEIT< matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        EIT
        elem_centers
        elem_volumes
        coeffs
        % zern_coeffs
        % cheb_coeffs
        % zern_set
        % cheb_set
        cond_values
        coeff_matrix
        saved_info_elem_currents
        n_coeffs
        img
        img_0
        elem_currents_0
        % elem_currents
        phantom
        unnormalized_coil_voltages
        pert_amplitude
        coil_system
        unnormalized_coil_voltages_0
    end

    properties(Hidden)
        scaled_elem_centers
        elem_currents_
        % for storage of elem_currents
        already_loaded_cycle
        already_loaded_elem_currents
    end

    methods
        function obj = InverseProblemEIT(phantom, current_ampl, ...
                coeffs)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % Avoid targets in this model
            if isfield(phantom, 'extra')
                phantom = rmfield(phantom,'extra');
            end

            % use unitary background
            phantom.background = 1;
            obj.phantom = phantom;
            
            obj.EIT = EIT(phantom, current_ampl);

            obj.img_0 = obj.EIT.img;
            obj.img = obj.EIT.img;

            obj.elem_centers = obj.EIT.elem_centers;
            obj.elem_volumes = obj.EIT.elem_volumes;

            obj.scaled_elem_centers = helpers.cylindrical_elem_centers(obj.elem_centers);
            
            obj.n_coeffs = obj.assign_coeffs(coeffs);
        end
    end
        
    methods (Abstract)
        n_coeffs = assign_coeffs(obj, coeffs)
        make_basis(obj)
        calc_cond_values(obj, pert_amplitude)
    end

    methods
        function assign_values(obj, coeff_to_see)
            coeff_to_see = coeff_to_see(:);

            idx = find(ismember(obj.coeff_matrix', coeff_to_see', 'row'));

            if isempty(idx)
                error('Wrong coefficients')
            end

            obj.assign_cond_values_from_idx(idx)
        end

        function assign_cond_values_from_idx(obj, cond_idx)
            if cond_idx ~= 0
                obj.img.elem_data = obj.cond_values(:,cond_idx);
                obj.EIT.img.elem_data = obj.cond_values(:,cond_idx);
            else
                obj.img.elem_data = obj.img_0.elem_data;
                obj.EIT.img.elem_data = obj.img_0.elem_data;
            end
        end

        function load_saved_info_elem_currents(obj, fullpath)
            if nargin == 1
                [file,location] = uigetfile('*.mat');
                selectedfile = fullfile(location,file);
                data = load(selectedfile);
                obj.saved_info_elem_currents = data;
            else
                data = load(fullpath);
                obj.saved_info_elem_currents = data;
            end
        end

        function elem_currents = get_elem_currents(obj, idx)
            arguments
                obj
                idx = [] % it is only for the cond_value
            end
            if isempty(obj.saved_info_elem_currents)
                if isempty(idx)
                    elem_currents = obj.elem_currents_;
                else
                    elem_currents = obj.elem_currents_(:,:,:,idx);
                end
            else
                % find in which file that index is
                for ii = 1:obj.saved_info_elem_currents(1).n_cycles
                    l = find(idx == obj.saved_info_elem_currents(ii).cond_values_idxs,1);
                    if ~isempty(l)
                        if ii == obj.already_loaded_cycle
                            % this cycle was loaded already
                            
                        else % we load the file first
                            fullpath = obj.saved_info_elem_currents(ii).dirname + '/' + obj.saved_info_elem_currents(ii).filenames;
                            buffer_file = matfile(fullpath);
                            obj.already_loaded_elem_currents = buffer_file.elem_currents;
                            obj.already_loaded_cycle = ii;
                        end
                        elem_currents = obj.already_loaded_elem_currents(:,:,:,l);
                        break
                    end
                end
            end
        end

        function save_all_currents(obj, dirname, cond_values_per_cycle)
            arguments
                obj
                dirname = string(datetime('now','format', 'yyyyMMdd_HH_mm_ss'))
                cond_values_per_cycle = 50
            end

            measurements_idx = 1:size(obj.EIT.volt_strct.volt,2); % 1:obj.phantom.n_elec
            % obj.elem_currents = zeros(length(obj.elem_centers), 3, obj.phantom.n_elec, size(obj.cond_values,2));
             
            n_cond_values = size(obj.cond_values,2);
            n_cycles = ceil(n_cond_values/cond_values_per_cycle);

            filenames = cell(n_cycles,1);
            cond_values_idxs = cell(n_cycles,1);
            l0=0;
            for ii = 1:n_cycles
                filenames{ii} = sprintf('%03d.mat', ii);
                idxs = l0+(1:cond_values_per_cycle);
                l0 = idxs(end);
                if l0 > n_cond_values
                    idxs = idxs(1): n_cond_values;
                end
                cond_values_idxs{ii} = idxs;             
            end
            mkdir(dirname);
            % info
            obj.saved_info_elem_currents = struct('dirname', dirname, ...
                'filenames', filenames, ...
                'cond_values_idxs', cond_values_idxs, ...
                'n_cond_values', n_cond_values, ...
                'n_cycles', n_cycles);

            fullpath = dirname + '/' + 'info_elem_currents.mat';
            save(fullpath, 'filenames', ...
                    'cond_values_idxs', ...
                    'n_cond_values', ...
                    'n_cycles', ...
                    'dirname')
            
            % unperturbated current
            obj.elem_currents_0 = obj.calc_elem_current(0,measurements_idx);

            % all perturbated currents
            progressbar(0,0)
            for ii = 1:n_cycles
                elem_currents = zeros(length(obj.elem_centers), 3, obj.phantom.n_elec, length(cond_values_idxs{ii})); %#ok<PROPLC>
                progressbar([],0) % Reset 2nd bar
                for ll = 1:length(cond_values_idxs{ii})
                    idx = cond_values_idxs{ii}(ll);
                    elem_currents(:,:,:,ll) = obj.calc_elem_current(idx,measurements_idx); %#ok<PROPLC>
                    progressbar([],ll/length(cond_values_idxs{ii})) % Update 2nd bar
                end
                fullpath = dirname + '/' + filenames{ii};
                save(fullpath, 'elem_currents', 'idx', '-v7.3');
                fprintf('%s Cycle %03d saved in %s\n', string(datetime('now','format', 'dd.MM.yy HH:mm:ss')), ii, fullpath)
                progressbar(ii/n_cycles) % Update 1st bar
            end
        end

        function calc_all_currents(obj)
            obj.elem_currents_ = zeros(length(obj.elem_centers), 3, obj.phantom.n_elec, size(obj.cond_values,2));
            measurements_idx = 1:size(obj.EIT.volt_strct.volt,2); % 1:obj.phantom.n_elec
            disp('Calculating currents for each perturbation...')
            tic
            progressbar
            for ii = 0:size(obj.cond_values,2)
                if ii == 0
                    obj.elem_currents_0 = obj.calc_elem_current(ii,measurements_idx);
                    continue
                end
                obj.elem_currents_(:,:,:,ii) = obj.calc_elem_current(ii,measurements_idx);
                waitbar(ii/size(obj.cond_values,2),f);   
                progressbar(ii/size(obj.cond_values,2))
            end
            t = duration(0,0,toc, 'Format', 'hh:mm:ss');
            fprintf('It took %s\n', t)
        end

        function elem_currents = calc_elem_current(obj, cond_idx, measurements_idx)
            assign_cond_values_from_idx(obj, cond_idx)
            obj.EIT.calc_elem_current(measurements_idx)
            elem_currents = obj.EIT.elem_currents;
        end

        function attach_coil_system(obj, coil_system)

            obj.coil_system = coil_system;

        end

        function [unnormalized_coil_voltages_0, unnormalized_coil_voltages] = calc_coil_integrals(obj)
            model = struct('elem_centers', [], ...
                           'elem_volumes', [], ...
                           'elem_curr', []);
            model.elem_centers = obj.elem_centers;
            model.elem_volumes = obj.elem_volumes;
            
            disp('Calculating integrals for each perturbation...')
            f = waitbar(0,'Calculating integrals for each perturbation...');
            unnormalized_coil_voltages = zeros(obj.EIT.n_elec^2,size(obj.cond_values,2));
            
            tic;
            for ii = 0:size(obj.cond_values,2)
                if ii == 0
                    model.elem_curr = obj.elem_currents_0;
                    unnormalized_coil_voltages_0 = reshape(obj.coil_system.calc_coil_integrals( ...
                                    model, 1:obj.EIT.n_elec), [], 1);
                    continue
                end
                model.elem_curr = obj.elem_currents(:,:,:,ii);
                unnormalized_coil_voltages(:,ii) = reshape(obj.coil_system.calc_coil_integrals( ...
                                    model, 1:obj.EIT.n_elec), [], 1);
                waitbar(ii/size(obj.cond_values,2),f);
            end
            close(f)
            t = duration(0,0,toc, 'Format', 'hh:mm:ss');
            fprintf('It took %s\n', t)

            obj.unnormalized_coil_voltages = unnormalized_coil_voltages;
            obj.unnormalized_coil_voltages_0 = unnormalized_coil_voltages_0;
        end

        function jacobian = calc_jacobian(obj)
            if isempty(obj.unnormalized_coil_voltages)
                obj.calc_coil_integrals();
            end

            d_unnormalized_voltages = obj.unnormalized_coil_voltages-obj.unnormalized_coil_voltages_0;
            jacobian = d_unnormalized_voltages/obj.pert_amplitude;

        end
    end
end