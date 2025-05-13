
function supp_dist_util_show_one_surface_map()

    params.SCRATCH_SPACE = 'E:\temp';
    
    [filename, pathname] = uigetfile(fullfile(params.SCRATCH_SPACE, '*.mat'));
    if isequal(filename,0), return; end
    
    load(fullfile(pathname, filename), 'config', 'params', ...
            'orig_V1', 'orig_V2', 'd_array', 'Ts', 'V1_reg', 'R_full', 'results');

    % -----------------------------------------------------------------------------------------
    % Show graphic of distortion mapped onto surface of moving volume...
    %
    params.SMOOTH_VOLUME_FOR_DISPLAY = false;
    dist_colour_map_on_volume_surface(config, params, d_array, Ts, V1_reg, R_full);
    %
    % -----------------------------------------------------------------------------------------

end