



function nifti_file_out = util_make_vp_nifti_vols(nifti_dir_in, dicom_dir_in)

    if nargin < 1
        
        nifti_dir_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001\Segmentations';
        dicom_dir_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001\Segmentations\7094_950_ADC_10_6_mm_s_Supine.nii.gz';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001\Segmentations\7094_1751_Water_DCE_Supine_LastPhase.nii.gz';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001\Segmentations\7094_2151_Water_DCE_Prone_FirstPhase.nii.gz';
        %nifti_file_in = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0001\Segmentations\7094_2050_ADC_10_6_mm_s_Prone.nii.gz';
        
    end
    
    nifti_dir_out = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\P0VP3\Segmentations';
    
    if ~exist(nifti_dir_out, 'dir'), mkdir(nifti_dir_out); end
    
    nii_files = dir(fullfile(nifti_dir_in, '*.nii.gz'));
    
    for n = 1:length(nii_files)
   
        nifti_file_in = nii_files(n).name;
       
        d = dir(fullfile(dicom_dir_in, '*.'));
        
        for i = 3:length(d)
            
            if contains(nifti_file_in, d(i).name)
                this_dicom_dir = fullfile(dicom_dir_in, d(i).name);
                break
            end
            
        end
        
        df = dir(fullfile(this_dicom_dir, '*.dcm'));
        
        H1 = dicominfo(fullfile(this_dicom_dir, df(1).name));
        Hn = dicominfo(fullfile(this_dicom_dir, df(end).name));
        
        pixel_spacing = H1.PixelSpacing;
        try
            slice_spacing = H1.SpacingBetweenSlices;
        catch
            slice_spacing = H1.SliceThickness;
        end
        
        if     H1.SliceLocation > Hn.SliceLocation
            slice_direction = 'IS';
        elseif H1.SliceLocation < Hn.SliceLocation
            slice_direction = 'SI';
        end
    
        nifti_file_out = fullfile(nifti_dir_out, ...
            strrep(nifti_file_in, '.nii.gz', '_VP.nii'));
        
        [VOLS, num_found_vols, H] = util_nifti_vol_to_masks(fullfile(nifti_dir_in, nifti_file_in));
        
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
            
            lo_slice = slices_present(1);
            hi_slice = slices_present(end);
            
            num_slices = hi_slice - lo_slice + 1;
            
            max_w = 20;
            max_h = 15;
            
            centre_slice = lo_slice + round( (hi_slice - lo_slice) / 2.0);
            
            centre_im = V(:,:,centre_slice);
            
            cen = regionprops(centre_im, 'Centroid').Centroid;
            
            if size(cen,2) > 1, cen = cen(1,:); end
            
            figure();
            
            count = 1; %0;
            
            switch slice_direction
                case 'SI'
                    slice_beg = lo_slice;
                    slice_end = hi_slice;
                    slice_inc = +1;
                case 'IS'
                    slice_beg = hi_slice;
                    slice_end = lo_slice;
                    slice_inc = -1;
            end
            
            for s = slice_beg:slice_inc:slice_end
                
                this_im = V(:,:,s);
                
                size_inc = (count / num_slices);
                
                w = (max_w * size_inc) / pixel_spacing(1);    %(5 + size_inc) * pixel_spacing(1);  % If not square pixels, check this...!
                h = (max_h * size_inc) / pixel_spacing(2);    %(1 + size_inc) * pixel_spacing(2);
              
                count = count + 1;
                 
                imshow(this_im);
                
                roi = drawellipse(gca, 'Center',cen, 'SemiAxes',[w,h]);
                
                mask = createMask(roi);
                
                imshow(mask);
                
                V(:,:,s) = mask;
            
                area(s) = numel(nonzeros(mask)) * pixel_spacing(1) * pixel_spacing(2); %#ok<AGROW>
                
                disp(['Pixels = ' num2str(numel(nonzeros(mask)))]);
                disp(['Area   = ' num2str(numel(nonzeros(mask)) * pixel_spacing(1) * pixel_spacing(2))]);
                
            end
            
            disp(' ');
            disp(['Slices = ' num2str(num_slices)]);
            disp(['Volume = ' num2str(sum(area) * slice_spacing)]); 
            disp(' ');
            
            if num_found_vols > 1
                
                VOL_OUT = VOL_OUT + data_fn_h(v * ones(size(V)) .* single(V));
                
            else
                
                VOL_OUT = data_fn_h(ones(size(V)) .* single(V));
                
            end
            
            clear area   % Important...!
            
        end
        
        VOL_OUT = permute(VOL_OUT, [2,1,3]);
        
        figure;
        for s = 1:size(VOL_OUT,3)
            imshow(VOL_OUT(:,:,s),[]);
            drawnow();
        end
        
        niftiwrite(VOL_OUT, nifti_file_out, H, 'Compressed',true);
        
    end
    
end
