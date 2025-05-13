


function [results, hf] = dist_calc_distortion(params, orig_V1, R1, orig_V2, R2)
% Apply the algorithm to quantify and map distortion...
%
% Input:
%
%   orig_V1 : 'moving volume' labelled here as the 'original' volume in space R1
%           : as a 3-D BW array (1 inside volume; 0 outside)
%   R1      : spatial reference object associated with V1 (e.g. as created by MATLAB function 'imref3d()')
%
%   orig_V2 : 'fixed volume'  in space R2
%           : as a 3-D BW array (1 inside volume; 0 outside)
%   R2      : spatial reference object associated with V2
%
%   params.SCRATCH_SPACE             : A path to a location for storing temporary files
%
%   params.PYRAMID_LEVELS            : number of pyramid levels to be used in registration algorithm
%   params.SLICES_BY_REGION          : E.g. [3,3,4, 1] number of slices in V2 to be considered as [APEX,MID,BASE, z-pad]
%                                    : 'z-pad' is the number of blank slices V2 is padded with at both extents
%   params.SMOOTH_VOLUME_FOR_DISPLAY : [true/false] whether to smooth for final map display (cosmetic effect only)
%   params.DISTORTION_MAP_MAX_MM     : the maximum value [mm] of distortion colour scale on map
%
% Output:
%
%   results : structure containing distortion results for this pair of volumes
%
%        hf : handle to figure showing the distortion colour-map
%

    % ---------------------------------------------------------------------------------------------------------------------------------
    % Configure...
    config.RADIAL_STEP       = 0.1;     % Resolution of radial step in search for boundaries (pixels)
    config.SEARCH_RADIUS_PIX = 100.0;   % Extent of radial search for boundaries (pixels)
    config.SURFACE_THRESHOLD = 0.5;     % Threshold for grey-level of boundary (images already scaled to have intensities 0 < I < 1)
    config.VERBOSE           = params.VERBOSE;%1;       % Produce verbose output or not
    % ---------------------------------------------------------------------------------------------------------------------------------

    % -------------------------------------------------------------------------------------------------------------------------------
    % Three step algorithm...
    
    % (1) Register volumes by means of a translation transform...
    %
    %      V1_reg is registered 'moving' volume in combined space
    %      R_full is the spatial reference object ('combined space') encompassing both volumes
    %      V2 changes: same volume but now interpolated to 'combined space'
    
    [V1_reg, V2, R_full, t_transl]   = dist_register_volumes(orig_V1, R1, orig_V2, R2, params.PYRAMID_LEVELS);

    % (2) Translate each slice so that slice centroids coincide...
    %
    %      V1 is the new 'moving' volume after slice centroid translations have been applied
    %      V1C holds the new slice centroid positions
    %      Ts holds the slice-by-slice translations required to map centroids together
    %
    [V1, V1C, Ts]                    = dist_translate_moving_volume_slices_by_centroid(V1_reg, V2, R_full, config.SURFACE_THRESHOLD);
    
    % (3) Compute local radial distortion vectors...
    %
    %      d and d_res hold the r.m.s. distortion and resultant distortion values in mm respectively
    %      d_array and d_array_res hold the distortion/resultant distortion values in mm by slice and radial angle 'phi'
    %
    [d, sem_d, d_res, sem_d_res, d_array, d_array_res] = dist_par_compute_radial_distortion_by_slice(config, V1, V1C, V2, R_full, Ts);
    %[d, sem_d, d_res, sem_d_res, d_array, d_array_res] = dist_compute_radial_distortion_by_slice(config, V1, V1C, V2, R_full, Ts);
 
    % -------------------------------------------------------------------------------------------------------------------------------
    
    
    % -----------------------------------------------------------------------------------------
    % Show graphic of distortion mapped onto surface of moving volume...
    %
    hf = dist_colour_map_on_volume_surface(config, params, d_array, Ts, V1_reg, R_full);
    %
    % -----------------------------------------------------------------------------------------

    
    
    % -----------------------------------------------------------------------------------------
    % Populate results structure...
    
    SR                   = params.SLICES_BY_REGION;
    
    pt_pos               = params.PT_POSITION;
    
    results.t_transl     = t_transl;
    results.t_rotate     = affine3d(eye(4));   % No rotation allowed at present - room for expansion of code
    
    results.d            = d;  
    results.sem_d        = sem_d;
    
    results.d_res        = d_res;
    results.sem_d_res    = sem_d_res;
    
    results.d_array      = d_array;
    results.d_array_res  = d_array_res;
    
    results.Ts           = Ts;
   
    results.d_base       = calc_rms_distortion(d_array, SR, 'BASE');
    results.d_mid        = calc_rms_distortion(d_array, SR, 'MID');
    results.d_apex       = calc_rms_distortion(d_array, SR, 'APEX');
    
    results.d_ant        = calc_rms_distortion(d_array, SR, 'ANTERIOR');
    results.d_post       = calc_rms_distortion(d_array, SR, 'POSTERIOR');
    
    results.d_base_res   = calc_rms_distortion(d_array_res, SR, 'BASE');
    results.d_mid_res    = calc_rms_distortion(d_array_res, SR, 'MID');
    results.d_apex_res   = calc_rms_distortion(d_array_res, SR, 'APEX');
    
    switch pt_pos
        
        % 'phi' goes clockwise from Q1 -> Q4 such that it starts posterior
        %  on supine imaging and anterior on prone imaging...
        
        case 'SUPINE'
            results.d_ant_res    = calc_rms_distortion(d_array_res, SR, 'Q3-Q4');
            results.d_post_res   = calc_rms_distortion(d_array_res, SR, 'Q1-Q2');
        case 'PRONE'
            results.d_ant_res    = calc_rms_distortion(d_array_res, SR, 'Q1-Q2');
            results.d_post_res   = calc_rms_distortion(d_array_res, SR, 'Q3-Q4');
        otherwise
            error('Not a valid patient position...!');
            
    end
    
    % -----------------------------------------------------------------------------------------

    % Optional line to save variables...
    save(fullfile(params.SCRATCH_SPACE, 'temp.mat'), 'config', 'params', ...
            'orig_V1', 'orig_V2', 'd_array', 'Ts', 'V1_reg', 'R_full', 'results');

end



function d_region = calc_rms_distortion(d_array, SR, region)
% Calculate r.m.s. distortion values by region (breaking down by slice
% position and radial angle)...

    pad_by_lo = SR(4);
    %pad_by_hi = SR(5);
    
    slices = pad_by_lo:pad_by_lo+SR(1)+SR(2)+SR(3);
    phi_range = 1:360;
    
    switch region
        case 'BASE' %%%'APEX'
            slices = pad_by_lo+1:pad_by_lo+SR(1);
        case 'MID'
            slices = pad_by_lo+SR(1)+1:pad_by_lo+SR(1)+SR(2);
        case 'APEX' %%%'BASE'
            slices = pad_by_lo+SR(1)+SR(2)+1:pad_by_lo+SR(1)+SR(2)+SR(3);
        case 'Q3-Q4'
            phi_range = 181:360;    
        case 'Q1-Q2'
            phi_range = 1:180;      
    end
    
    c = 1;
    for s = slices
        for phi = phi_range
            d(c) = d_array(phi, s); %#ok<AGROW>
            c = c + 1;
        end
    end
    
    d_region = my_rms(d(~isnan(d)));
    
end

