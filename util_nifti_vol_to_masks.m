


function [masks, num_found_masks, h] = util_nifti_vol_to_masks(nii_file_in, roi_numbers_in_nifti)

    V = niftiread(nii_file_in);
    
    h = niftiinfo(nii_file_in);
    
    dims = size(V);
    
    % Swap x- and y-axes: seems to be necessary empirically...
    % [ May be related to ROW - COL swap in matrix when read in from NIFTI file..? ] 
    
    if     numel(dims) == 3
        V = permute(V, [2,1,3]); 
    elseif numel(dims) == 4
        V = permute(V, [2,1,3,4]);
    else
        error('Dimensionality error in NIFTI volume!');
    end

    num_found_masks = 0;
    
    for r = 1:size(roi_numbers_in_nifti,2)
        
        roi_number = roi_numbers_in_nifti(r);

        this_mask = (V == roi_number);
        
        if numel(nonzeros(this_mask)) == 0
            
            warning(['Null mask found for ROI # : ' num2str(roi_number) ]);
            %error(['Didn''t find a valid mask for ROI # : ' num2str(roi_number) ]);
            
        %else
        
        end
           
        num_found_masks = num_found_masks + 1;
            
        masks(:,:,:,num_found_masks) = this_mask; %#ok<AGROW>
         
        %end
        
    end
    
    
    %{
    found_roi = true;
    
    roi_number = 0;
    
    num_found_masks = 0;
    
    while found_roi
        
        if num_found_masks + 1 > size(roi_numbers_in_nifti,2), break; end
        
        roi_number = roi_number + 1;
        
        this_mask = (V == roi_number);
        
        if numel(nonzeros(this_mask)) == 0
            
            found_roi = false;
            
        else
            
            masks(:,:,:, roi_number) = this_mask;
    
            num_found_masks = roi_number;
            
        end
        
    end
    %}
    
    [~, f_name, ~] = fileparts(nii_file_in);
    
    disp(['In ' f_name ' found ' num2str(num_found_masks) ' ROI masks.']);
    
end