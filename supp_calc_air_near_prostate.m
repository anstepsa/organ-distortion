


function vol_air_cm3 = supp_calc_air_near_prostate(params, V, MASKS, voxel_size)

    % Calculate the volume of air within e.g., 10 mm of prostate...
    
    % Do this by dilating the prostate VOI by 10 mm and forming intersection 
    % of this with air VOI...

    vol_air_cm3 = -1;

    dims = size(MASKS);
    
    if size(dims) < 4
        warning('Only 1 ROI: can''t calculate air volume!');
        return
    end
    
    PROSTATE_MASK = squeeze(MASKS(:,:,:,1));
    AIR_MASK      = squeeze(MASKS(:,:,:,2));
    
    pixel_h = voxel_size(2);
    pixel_w = voxel_size(1);
    
    if pixel_h ~= pixel_w
        warning('Not square pixels: can''t dilate mask easily!');
        return
    end
    
    n = numel(params.DILATE_BY_MM);
    
    for i = 1:n
    
        disp(['Calculating air within ' num2str(params.DILATE_BY_MM(i)) ' mm of prostate!']);
        
        r = round(params.DILATE_BY_MM(i) / pixel_w);  % In pixels...

        vox_vol_mm3 = (pixel_h ^ 2) * voxel_size(3);

        tic
        
        dilated_prostate_mask = imdilate(PROSTATE_MASK, strel('sphere', r));%disk', r, 8));
        
        toc
        
        vol_air_mask = AIR_MASK .* dilated_prostate_mask;

        vol_air_cm3(:,i) = numel(nonzeros(vol_air_mask)) * vox_vol_mm3 / 1000;

        if params.VERBOSE

            % Show the air and original PROSTATE ROIs...

            hf = figure();
            set(hf, 'Color','w', 'Units','normalized', 'Position',[0.2,0.05,0.6,0.8]);

            for s = 1:size(V,3)
                I = double(V(:,:,s));
                R = double(AIR_MASK(:,:,s));
                A = double(vol_air_mask(:,:,s));
                P = double(PROSTATE_MASK(:,:,s));
                D = double(dilated_prostate_mask(:,:,s));
                im = I + (500 * P) + (750 * D) + (1000 * R) + (1250 * A);
                imshow(im, [], 'InitialMagnification',500);
                drawnow();
                if numel(nonzeros(A)) > 0
                    pause();
                end
            end

        end

    end
    
end