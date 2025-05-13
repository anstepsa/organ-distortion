


function T_RESULTS_BOTH = supp_dist_control_script_v4(config)

    % -- Set-up --------------------------------------------------------------
    
    config.MARGIN_PIXELS     = 50;
    
    config.PT_POSITIONS      = {'SUPINE', 'PRONE'};
    
    config.IMAGE_TYPES       = {'DIS', 'REF'};
     
    config.PT_BASE_DIR       = fullfile(config.BASE_DIR, config.PT_ID);
   
    config.SEGMENTATION_DIR  = fullfile(config.PT_BASE_DIR, 'Segmentations');
    
    % ------------------------------------------------------------------------
    
    switch config.PT_ID
        case  'V0000'
            n_pg = 1; n_rectum = 1; n_air = 1;  
        case  'V0001'   % [ was V0002 in previous data run ]
            n_pg = 1; n_rectum = 3; n_air = 2; 
        otherwise   % i.e., the default...
            n_pg = 1; n_rectum = 2; n_air = 3;  
    end
        
                                  %  Voxel value in NIFTI for associated ROIs
                                  %    PROSTATE, [PROSTATE,AIR,RECTUM]
    config.ROI_NUMBERS_IN_NIFTIS = { {    1,     [n_pg,    n_air,    n_rectum]}; ...   % {SUPINE: DIS, REF}
                                     {    1,     [n_pg,    n_air,    n_rectum]} };     % {PRONE : DIS, REF}
                                  
                                 % [ Could be made pt. dependant... ]
                                  
    config = get_patient_info(config);
     
    % -- End set-up ----------------------------------------------------------  
    
    % Open and begin results files...
    
    p_results_dir = fullfile(config.OUT_DIR, config.PT_ID); 
    if ~exist(p_results_dir, 'dir'), mkdir(p_results_dir); end
    
    T_RESULTS_BOTH = [];
    
    % Loop over SUPINE and PRONE cases...
    
    for pt_position_index = 1:2%:1
   
        patient_position = config.PT_POSITIONS{pt_position_index};
           
        roi_numbers_in_niftis = config.ROI_NUMBERS_IN_NIFTIS{pt_position_index};
        
        % Loop over DIS-torted and REF-erence images...
        
        for im_type_index = 1:2
            
            im_type = config.IMAGE_TYPES{im_type_index};
            
            % Read in DICOM files into volume(s) and ROI volume from NIFTI...
            
            dicom_dir.(im_type) = fullfile(config.PT_BASE_DIR, config.([patient_position '_DICOM_DIRS']){im_type_index});
        
            switch im_type
                case 'DIS'
                    roi_files.(im_type) = { fullfile(config.SEGMENTATION_DIR, config.([patient_position '_ROI_FILES']){1}) };
                case 'REF'
                    roi_files.(im_type) = { fullfile(config.SEGMENTATION_DIR, config.([patient_position '_ROI_FILES']){2}), ...
                                            fullfile(config.SEGMENTATION_DIR, config.([patient_position '_ROI_FILES']){3}), ...
                                            fullfile(config.SEGMENTATION_DIR, config.([patient_position '_ROI_FILES']){4}) };
            end
            
            [dicom_vol.(im_type),  hdr.(im_type), slice_position.(im_type)] ...
                                                   = util_read_dicom_dir(dicom_dir.(im_type));
            
            roi_nums = roi_numbers_in_niftis{im_type_index};
            
            [MASKS.(im_type), num_masks.(im_type)] = supp_read_masks_from_niftis(im_type, roi_files.(im_type), roi_nums);
           
                            % This always returns the PROSTATE mask in 4th index = 1
                            %       "             the AIR      mask in 4th index = 2
                         
            
            % Get positional information from DICOM header...
            
            stats.PIXEL_SPACING   = hdr.(im_type){1}.PixelSpacing;
            stats.IMG_ORIENTATION = hdr.(im_type){1}.ImageOrientationPatient;
            
            try
                stats.SLICE_SPACING = hdr.(im_type){1}.SpacingBetweenSlices;
            catch
                stats.SLICE_SPACING = hdr.(im_type){1}.SliceThickness;
            end
            
            
            % If first slice is pos < last slice : we're imaging IS  <- reverse the MASK, slice and DICOM information...
            % If first_slice is pos > last_slice : we're imaging SI  <- normal case
            
            if (slice_position.(im_type)(1,3) < slice_position.(im_type)(end,3)) 
                [MASKS.(im_type), dicom_vol.(im_type), slice_position.(im_type)] = ...
                    util_reverse_slice_order_in_vols(MASKS.(im_type), dicom_vol.(im_type), slice_position.(im_type));           
            end
            
            
            sp.(im_type)    = [stats.PIXEL_SPACING(1), stats.PIXEL_SPACING(2), stats.SLICE_SPACING];
                
            R.(im_type)     = util_construct_spatial_ref_object(MASKS.(im_type)(:,:,:,1), stats, ...
                                                   sp.(im_type), slice_position.(im_type));
                                               
                                             % Will always be at least one MASK, so use '1' here.
                                             % Spatial reference object is the same for all masks in
                                             % this image type...
            
            % Render pictures and form variable 'V'...
            for m = 1:num_masks.(im_type)
                
                hf = figure; set(hf, 'Color','w');
                
                mask = MASKS.(im_type)(:,:,:,m);
                
                % Show the volume in a figure in 3-D...
                supp_dist_visualise_volume(false, false, mask, [stats.PIXEL_SPACING(1), stats.PIXEL_SPACING(2), stats.SLICE_SPACING]); 
                
                savefig(hf, fullfile(p_results_dir, ['Volume_Render_' patient_position '_' im_type '_' num2str(m) '.fig']));
               
                if config.VERBOSE
                    % Sanity check : show slice by slice in 2-D...
                    supp_show_mask_and_dicom(mask, dicom_vol.(im_type));
                end
                
                V.(im_type){m}  = single(mask);
                
                % Calculate area of mask by slice...
                mask_area = [];
                for ss = 1:size(mask, 3)
                    mask_area(ss,1) = numel(nonzeros(mask(:,:,ss))); %#ok<AGROW>
                end
                mask_area = nonzeros(mask_area);
                mask_area = mask_area * stats.PIXEL_SPACING(1) * stats.PIXEL_SPACING(2);  
                
                % Crop the volumes for ease and speed of matching: iterate over
                % suitable number of margin pixels (in 'x' and 'y')...
                
                if ~isempty(mask_area) && any(mask_area > 0)
                    for margin = config.MARGIN_PIXELS:-5:0
                        [Vc, Rc, ok] = supp_dist_crop_volume(V.(im_type){m}, R.(im_type), margin);
                        if ok > 0
                            V.(im_type){m} = Vc; RO.(im_type){m} = Rc;
                            break
                        end
                    end
                end
                
                % Calculate volume of organ from areas...
                total_volume    = supp_dist_simpsons_rule_area([0.0 mask_area' 0.0], stats.SLICE_SPACING);
                total_ni_volume = sum(mask_area) * stats.SLICE_SPACING;
          
                results.(im_type).INTERP_VOLUME(m) = total_volume    / 1000.0;
                results.(im_type).ROUGH_VOLUME(m)  = total_ni_volume / 1000.0;
                
            end
            
            % If we have PG in '1' and 'AIR' in '2' then show this plot too...
            if num_masks.(im_type) == 2
            
                hf = figure; set(hf, 'Color','w');

                mask = MASKS.(im_type)(:,:,:,1) + MASKS.(im_type)(:,:,:,2);

                % Show the volume in a figure in 3-D...
                supp_dist_visualise_volume(false, false, mask, [stats.PIXEL_SPACING(1), stats.PIXEL_SPACING(2), stats.SLICE_SPACING]); 
                savefig(hf, fullfile(p_results_dir, ['Volume_Render_' patient_position '_' im_type '_1_and_2.fig']));
                
            end
            
            % Calculate air volume intersecting dilated PG...
            
            params.DILATE_BY_MM = [10, 15, 20];
            params.VERBOSE      = config.VERBOSE;
            
            if config.DO_CALC_AIR_VOLUME
            
                air_volume_cm3 = supp_calc_air_near_prostate(params, dicom_vol.(im_type), MASKS.(im_type), ...
                    [stats.PIXEL_SPACING(1), stats.PIXEL_SPACING(2), stats.SLICE_SPACING]); 
            
            else
            
                air_volume_cm3 = -1 * ones(size(params.DILATE_BY_MM));
                
            end
            
            results.(im_type).AIR_VOLUME = air_volume_cm3;
                
        end
              
        % ------------------------------------------------------------------------------------------
        if config.DO_CALC_DISTORTION
       
            params.SCRATCH_SPACE             = config.WORK_DIR; 
            
            %: A path to a location for storing temporary files
            
            params.PYRAMID_LEVELS            = 2;     %: number of pyramid levels to be used in registration algorithm
            
            params.SLICES_BY_REGION          = get_volume_slices(V.REF{1});
            
            %: E.g. [3,3,4, pad_z_lo, pad_z_hi] number of slices in V2 to be considered as [APEX,MID,BASE, z-pad]
            %: 'z-pad' is the number of blank slices V2 is padded with at both extents
            
            params.SMOOTH_VOLUME_FOR_DISPLAY = false; %: [true | false] whether to smooth for final map display (cosmetic effect only)
            
            params.DISTORTION_MAP_MAX_MM     = 15;    %: the maximum value [mm] of distortion colour scale on map
            
            params.VERBOSE                   = config.VERBOSE;
            
            params.PT_POSITION               = patient_position;
            
            % -- Calculate distortion measures --------------------------------------------------------
            
            disp('Calculating distortion measures...');
            
            [d_stats, hf_out] = dist_calc_distortion(params, V.DIS{1}, RO.DIS{1}, V.REF{1}, RO.REF{1});
            
            % -----------------------------------------------------------------------------------------

            % Accrue distortion results...        
            results.DISTORTION = d_stats;
            
            save_str = patient_position;
            savefig(hf_out, fullfile(p_results_dir, ['Volume_Surface_Dist_Map_' save_str '.fig']));
            close(hf_out);
      
            mat_file = fullfile(params.SCRATCH_SPACE, 'temp.mat');
            movefile(mat_file, fullfile(p_results_dir, ['Proj_Dist_Vars_' save_str '.mat']));

        else
            
            results = fill_stats_with_nans(results);
    
        end
       
        T_RESULTS_BOTH = add_line_to_table_out(patient_position, config.PT_ID, T_RESULTS_BOTH, results);
        
    end
    
    return
   
    
    
    function T_temp = set_up_table_out(pos, p_id, stats)
        
        % ! Could save more output here ... !
        
        if isfield(stats.DISTORTION, 't_rotate')
            [ax, ang] = get_angle_and_axis(stats.DISTORTION.t_rotate.T);
        else
            ax = NaN(1,3); ang = NaN;
        end
        
        stats.DISTORTION.ax  = ax;
        stats.DISTORTION.ang = ang;
        
        if isfield(stats.DISTORTION, 't_transl')
            stats.DISTORTION.t_x = stats.DISTORTION.t_transl.T(4,1);
            stats.DISTORTION.t_y = stats.DISTORTION.t_transl.T(4,2);
            stats.DISTORTION.t_z = stats.DISTORTION.t_transl.T(4,3);
        else
            stats.DISTORTION.t_x = NaN;
            stats.DISTORTION.t_y = NaN;
            stats.DISTORTION.t_z = NaN;
        end
        
        % This sets up a one-row (temporary) table with the relevant output in it...
        T_temp = table(string(p_id), string(pos), ...
                          stats.DIS.INTERP_VOLUME, ...
                          stats.REF.INTERP_VOLUME(1), ...
                          stats.REF.INTERP_VOLUME(2), ...
                          stats.REF.AIR_VOLUME(1), ...
                          stats.REF.AIR_VOLUME(2), ...
                          stats.REF.AIR_VOLUME(3), ...
                          stats.DISTORTION.ax(1), ...
                          stats.DISTORTION.ax(2), ...
                          stats.DISTORTION.ax(3), ...
                          stats.DISTORTION.ang, ... 
                          stats.DISTORTION.t_x, ...
                          stats.DISTORTION.t_y, ...
                          stats.DISTORTION.t_z, ...                                              
                          stats.DISTORTION.d, ...
                          stats.DISTORTION.sem_d, ...
                          stats.DISTORTION.d_res, ...
                          stats.DISTORTION.sem_d_res, ...
                          stats.DISTORTION.d_base, ...
                          stats.DISTORTION.d_mid, ...
                          stats.DISTORTION.d_apex, ...
                          stats.DISTORTION.d_ant, ...
                          stats.DISTORTION.d_post, ...
                          stats.DISTORTION.d_base_res, ...
                          stats.DISTORTION.d_mid_res, ...
                          stats.DISTORTION.d_apex_res, ...
                          stats.DISTORTION.d_ant_res, ...
                          stats.DISTORTION.d_post_res, ...
                     'VariableNames',{'pt_id', 'scan_position', ...
                                      'adc_prostate_volume', ...
                                      'ref_prostate_volume', ...
                                      'ref_air_volume', ...
                                      'ref_prox_air_volume_1', ...
                                      'ref_prox_air_volume_2', ...
                                      'ref_prox_air_volume_3', ...
                                      'rot_ax_x',   'rot_ax_y', 'rot_ax_z', ...
                                      'rot_degrees', ...
                                      'trans_x',    'trans_y',  'trans_z', ...
                                      'd',          'sem_d', ...
                                      'd_res',      'sem_d_res', ...
                                      'd_base',     'd_mid',    'd_apex', ...
                                      'd_ant',      'd_post', ...
                                      'd_base_res', 'd_mid_res','d_apex_res', ...
                                      'd_ant_res',  'd_post_res', ...
                                      } ...
                       );
 
                   
                   %NaN, ...
                   %NaN, ...
                   
                   %stats.REF.AIR_VOLUME(2), ...
                   %stats.REF.AIR_VOLUME(3), ...
                   
    end
    
   
    function stats = fill_stats_with_nans(stats)
        
        stats.DISTORTION.ax(1)      = NaN;
        stats.DISTORTION.ax(2)      = NaN;
        stats.DISTORTION.ax(3)      = NaN;
        stats.DISTORTION.ang        = NaN;
        stats.DISTORTION.t_x        = NaN;
        stats.DISTORTION.t_y        = NaN;
        stats.DISTORTION.t_z        = NaN;
        stats.DISTORTION.d          = NaN;
        stats.DISTORTION.sem_d      = NaN;
        stats.DISTORTION.d_res      = NaN;
        stats.DISTORTION.sem_d_res  = NaN;
        stats.DISTORTION.d_base     = NaN;
        stats.DISTORTION.d_mid      = NaN;
        stats.DISTORTION.d_apex     = NaN;
        stats.DISTORTION.d_ant      = NaN;
        stats.DISTORTION.d_post     = NaN;
        stats.DISTORTION.d_base_res = NaN;
        stats.DISTORTION.d_mid_res  = NaN;
        stats.DISTORTION.d_apex_res = NaN;
        stats.DISTORTION.d_ant_res  = NaN;
        stats.DISTORTION.d_post_res = NaN;
            
    end


    % Given a transform matrix, extract a rotational axis and angle...
    function [ax, ang] = get_angle_and_axis(M)

        M = M(1:3,1:3);

        lambda = eig(M);
        [v, ~] = eig(M);

        col = find( (abs( real(lambda) - 1) < 0.0001), 1, 'first');
        ax = v(:,col)';

        a = trace(M);
        ang = acosd( (a-1) / 2 );

    end


    function T_OUT = add_line_to_table_out(pos, p_id, T_OUT, stats)
       
        % Start the output table if blank or add the new line to it...
        
        T_temp = set_up_table_out(pos, p_id, stats);
        
        if isempty(T_OUT)
            T_OUT = T_temp;
        else
            T_OUT = cat(1, T_OUT, T_temp);
        end
        
    end

      
 

    function [V, hdr, slice_pos] = util_read_dicom_dir(dicom_dir)

        [~, dir_name, ~] = fileparts(dicom_dir);
        
        df = dir(fullfile(dicom_dir, '*.dcm'));

        ndf = numel(df);

        try
            hdr{1} = dicominfo(fullfile(dicom_dir, df(1).name));
        catch
            error(['No DICOM files found in ' dicom_dir ' !']);
        end
        
        disp(['Reading ' num2str(ndf) ' DICOM files in ' dir_name]);
      
        for s = 1:ndf
            
            hdr{s} = dicominfo(fullfile(dicom_dir, df(s).name)); %#ok<AGROW>
            V(:,:,s) = dicomread(hdr{s}); %#ok<AGROW>

            slice_pos(:,s) = hdr{s}.ImagePositionPatient; %#ok<AGROW>
          
        end

        slice_pos = slice_pos';
        
    end


    function [v, dv, sl] = util_reverse_slice_order_in_vols(volume_matrix, dicom_vol, slice_position) 
        
        % If slices acquired in opposite direction, swap the order...
        % z must go 'Increasing' from 'Inferior -> Superior'... CHECK!!! ok
        % I.e. Volume displayed and held as prostate base ('superior') at high 'z'

        disp('Swapping z-order for slices!');
        
        dims = size(volume_matrix);
        
        if numel(dims) > 3
            num_vols = size(volume_matrix, 4);
        else
            num_vols = 1;
        end
        
        num_slices = size(volume_matrix, 3);
        
        v = zeros(size(volume_matrix));
        
        for i = 1:num_vols
        
            vol = volume_matrix(:,:,:,i);
            
            sl = zeros(num_slices,3);
          
            for s = 1:num_slices
                
                v(:,:,s,i) = vol(:,:,num_slices-s+1); 
                
                dv(:,:,s)  = dicom_vol(:,:,num_slices-s+1); %#ok<AGROW>
                
                sl(s,:) = slice_position(num_slices-s+1,:); 
                
            end
            
            % clear dv...? NO, not necessary...
            
        end
        
        v = squeeze(v);
        
    end

    function R = util_construct_spatial_ref_object(v, stats, sp, slice_position)
    
        % Construct a spatial reference object: see Matlab help on registration
        % techniques... (This relates the volume matrix to real-world position and
        % coordinate system.)
     
        slice_position = slice_position';
        
        io = stats.IMG_ORIENTATION;

        XWorldLimits = [slice_position(1,1), slice_position(1,1) + size(v,1) * sp(1)] + (sp(1) / 2.0);
        YWorldLimits = [slice_position(2,1), slice_position(2,1) + size(v,2) * sp(2)] + (sp(2) / 2.0);

        %{
        if (slice_position(3,end) - slice_position(3,1)) < 0 %> 0
            warning('Slices still listed in decreasing ''z''!');
        end
        %}
        
        ZWorldLimits = [slice_position(3,1), slice_position(3,1) + size(v,3) * sp(3)] + (sp(3) / 2.0);

        if ~( myequal(abs(io(1)), 1, 1e-2) && myequal(abs(io(5)), 1, 1e-2) )
            warning('Volume may have been acquired in oblique plane!');
        end

        R = imref3d(size(v), XWorldLimits, YWorldLimits, ZWorldLimits);
        
    end


    function tf = myequal(a, b, tol)

        tf = (abs( (a - b) / a ) < tol);

    end


    function SLICES_BY_REGION = get_slices_by_region(total_organ_slices)

        % Split available slices into 3 regions...
        
        nr = floor(total_organ_slices / 3);

        switch mod(total_organ_slices, 3)
            case 0
                SLICES_BY_REGION = [nr,nr,nr];
            case 1
                SLICES_BY_REGION = [nr,nr+1,nr];
            case 2
                SLICES_BY_REGION = [nr,nr+1,nr+1];
        end

    end

    function SLICES_BY_REGION = get_volume_slices(V)
      
        % Split volume into border (top, bottom) and three prostate regions...
        
        organ_slices = [];
        
        for s = 1:size(V,3)
            im = V(:,:,s);
            if any(im(:))
                organ_slices = cat(2, organ_slices, 1);
            else
                organ_slices = cat(2, organ_slices, 0);
            end
        end
        
        pad_z_lo =             find(organ_slices, 1, 'first') - 1;
        pad_z_hi = size(V,3) - find(organ_slices, 1, 'last');
        
        total_slices = sum(organ_slices);
        
        SLICES_BY_REGION = get_slices_by_region(total_slices);
        
        SLICES_BY_REGION = cat(2, SLICES_BY_REGION, [pad_z_lo, pad_z_hi]);
        
    end


    function config = get_patient_info(config)
        
        sd = dir(fullfile(config.PT_BASE_DIR, '*.'));
    
        nsd = length(sd);
        c = 1;
        for i = 1:nsd
            if ~contains(sd(i).name, '.') && ~contains(sd(i).name, 'Segmentations')
                these_dicom_dirs{c} = sd(i).name;   %#ok<AGROW>
                c = c + 1;
            end
        end
        
        config.SUPINE_DICOM_DIRS = get_match(these_dicom_dirs, 'SUPINE', {{'ADC'},{'DCE'}});
        config.PRONE_DICOM_DIRS  = get_match(these_dicom_dirs, 'PRONE',  {{'ADC'},{'DCE'}});
        
        rf = dir(fullfile(config.SEGMENTATION_DIR, '*.nii.gz'));
        
        nrf = length(rf);
        
        assert(nrf == 8);
        
        for i = 1:nrf
            these_roi_files{i} = rf(i).name; %#ok<AGROW>
        end
        
        % A {string} or {set of strings} that defines : i.e.: ADC-PROSTATE,  WATER-PROSTATE,  WATER-AIR,  WATER-RECTUM
        config.SUPINE_ROI_FILES = get_match(these_roi_files, 'SUPINE', ...
                                          {{'ADC','PG'}, {'DCE','PG'}, {'DCE','AIR'}, {'DCE','RECTUM'}});
        config.PRONE_ROI_FILES  = get_match(these_roi_files, 'PRONE',  ...
                                          {{'ADC','PG'}, {'DCE','PG'}, {'DCE','AIR'}, {'DCE','RECTUM'}});
        
        return
        
        
        function str_out = get_match(C, pos, im_type)
            str_out = cell([numel(im_type),1]);
            for ii = 1:length(im_type)
                num_matches = 0;
                for cc = 1:length(C)
                    if contains(C{cc}, pos, 'IgnoreCase',true)
                        found_str = zeros(numel(im_type{ii}),1);
                        for jj = 1:numel(im_type{ii})
                            if ~contains(C{cc}, im_type{ii}{jj}, 'IgnoreCase',true)
                                break
                            else
                                found_str(jj) = 1;
                            end
                        end
                        if all(found_str)
                            match_vals{ii}{num_matches+1} = C{cc}; %#ok<AGROW>
                            num_matches = num_matches + 1;
                        end
                    end
                end
                if num_matches < 1
                    warning('No matches in search for files or directories!');
                    str_out{ii} = '<no matches>';               
                elseif num_matches > 1
                    warning('More than one match in search for files or directories!');
                    str_out{ii} = '<not resolved unambiguously>';
                else
                    str_out{ii} = match_vals{ii}{1};
                end
            end
        end
             
    end

end

   %{
    % Print a string to screen and file...
    function output_line(fid, s)
        fprintf(1,   '%s\n', s);
        fprintf(fid, '%s\n', s);
    end
    %}

   %{
    config.SUPINE_DICOM_DIRS = {'950_ADC_10_6_mm_s_Supine', ...
                                '1751_Water_DCE_Supine_LastPhase'};
    
    config.PRONE_DICOM_DIRS  = {'2050_ADC_10_6_mm_s_Prone', ...
                                '2151_Water_DCE_Prone_FirstPhase'};
    
    config.SUPINE_ROI_FILES  = {'7094_950_ADC_10_6_mm_s_Supine.nii.gz', ...
                                '7094_1751_Water_DCE_Supine_LastPhase.nii.gz'};
    
    config.PRONE_ROI_FILES   = {'7094_2050_ADC_10_6_mm_s_Prone.nii.gz', ...
                                '7094_2151_Water_DCE_Prone_FirstPhase.nii.gz'};
    
    switch config.PT_ID
        case {'P0VP1', 'P0VP2', 'P0VP3'}
            temp_str = '_VP';
        case 'P0000'
            temp_str = '_marked';
        otherwise
            temp_str = '';
    end
    
    config.SUPINE_ROI_FILES = cellfun(@strrep, config.SUPINE_ROI_FILES, repmat({'.nii.gz'},[1,2]), repmat({[temp_str '.nii.gz']},[1,2]), 'UniformOutput',false);
    config.PRONE_ROI_FILES  = cellfun(@strrep, config.PRONE_ROI_FILES,  repmat({'.nii.gz'},[1,2]), repmat({[temp_str '.nii.gz']},[1,2]), 'UniformOutput',false);
    %}