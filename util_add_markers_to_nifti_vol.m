



function nifti_file_out = util_add_markers_to_nifti_vol(nifti_file_in)

    if nargin < 1
        
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0000\Segmentations\7094_950_ADC_10_6_mm_s_Supine.nii.gz';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0000\Segmentations\7094_1751_Water_DCE_Supine_LastPhase.nii.gz';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0000\Segmentations\7094_2151_Water_DCE_Prone_FirstPhase.nii.gz';
        nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0000\Segmentations\7094_2050_ADC_10_6_mm_s_Prone.nii.gz';
        
    end
    
    nifti_file_out = strrep(nifti_file_in, '.nii.gz', '_marked.nii');

    [VOLS, num_found_vols, H] = util_nifti_vol_to_masks(nifti_file_in);
  
    data_type = H.Datatype;
    data_fn_h = str2func(data_type);
      
    v_size = size(VOLS);
    
    VOL_OUT = data_fn_h(zeros(v_size(1:3)));
    
    for v = 1:num_found_vols

        if num_found_vols > 1
            V = VOLS(:,:,:,v);
        else
            V = VOLS;
        end
        
        slices_present = [];
        
        for s = 1:size(V,3)
            
            this_im = V(:,:,s);
            
            if any(this_im(:))
                slices_present = cat(1, slices_present, s);
            end
                       
        end
        
        top_im = V(:,:,slices_present(1));
        low_im = V(:,:,slices_present(end));
        
        top_cen = regionprops(top_im, 'Centroid').Centroid;
        low_cen = regionprops(low_im, 'Centroid').Centroid;
        
        if size(top_cen,2) > 1, top_cen = top_cen(1,:); end
        if size(low_cen,1) > 1, low_cen = low_cen(1,:); end
        
        top_R = round(top_cen(1));  low_R = round(low_cen(1));
        top_C = round(top_cen(2));  low_C = round(low_cen(2));
  
        width = size(V,1);
       
        top_slice = slices_present(1)   - 1;
        low_slice = slices_present(end) + 1;
        
        if top_slice < 1,         top_slice = 1;         end
        if low_slice > size(V,3), low_slice = size(V,3); end
        
        width_hi = round(width / 10);
        width_lo = round(width / 20);
        
        V(  (top_C - width_lo) : (top_C + width_lo), ...
            (top_R - width_hi) : (top_R + width_hi), ...
             top_slice) = true;
         
        V(  (low_C - width_hi) : (low_C + width_hi), ...
            (low_R - width_lo) : (low_R + width_lo), ...
             low_slice) = true;
       
        if num_found_vols > 1
            
            VOL_OUT = VOL_OUT + data_fn_h(v * ones(size(V)) .* single(V));
            
        else
            
            VOL_OUT = data_fn_h(ones(size(V)) .* single(V));
            
        end
       
    end
     
    VOL_OUT = permute(VOL_OUT, [2,1,3]);
    
    figure;
    for s = 1:size(VOL_OUT,3)
        imshow(VOL_OUT(:,:,s),[]);
        drawnow();
    end
    
    niftiwrite(VOL_OUT, nifti_file_out, H, 'Compressed',true);
    
end
