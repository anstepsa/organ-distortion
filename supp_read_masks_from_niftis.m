

function [masks, num_masks] = supp_read_masks_from_niftis(im_type, roi_niftis, roi_num_in_nifti)

    nr = length(roi_niftis);
    
    switch im_type
        case 'DIS'
            assert(nr == 1);   % One ROI file only: 
                               %     1 : DIS : PROSTATE
        case 'REF'
            assert(nr == 3);   % Three ROI files:
                               %     1 : REF : PROSTATE
                               %     2 : REF : RECTUM
                               %     3 : REF : AIR
        otherwise
            error('Unknown ''im_type''...!');
    end
    
    for r = 1:nr
        
        [mask(:,:,:,r), num_found_masks] = util_nifti_vol_to_masks(roi_niftis{r}, roi_num_in_nifti(r));  %#ok<AGROW> % Must be only one ROI in file and we posit that it has value '1'...
          
        assert(num_found_masks == 1);
        
    end

    masks(:,:,:,1) = squeeze(mask(:,:,:,1));
    
    if     nr == 1
     
        num_masks = 1;
    
    elseif nr == 3
        
        AIR_MASK = single(squeeze(mask(:,:,:,2)));
        REC_MASK = single(squeeze(mask(:,:,:,3)));
        
        REC_AIR_MASK = AIR_MASK .* REC_MASK;
    
        masks(:,:,:,2) = logical(REC_AIR_MASK);
        
        num_masks = 2;
        
    end
    
end