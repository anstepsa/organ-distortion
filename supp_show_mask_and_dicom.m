


function supp_show_mask_and_dicom(mask, V)
            
    ts = datestr(datetime('now'), 'yyyy-mm-dd-HH-MM-SS'); %#ok<NASGU,DATST>
            
    num_slices = size(V, 3);
    
    hf = figure();
    set(hf, 'Color','w', 'Units','normalized', 'Position',[0.05,0.05,0.9,0.9]);
    
    for s = 1:num_slices
        
        subplot(1,2,1);
        imshow(V(:,:,s),[]);
      
        BW = mask(:,:,s);
        roi = boundarymask(BW, 4); %#ok<NASGU>
       
        subplot(1,2,2);
        
        I = double(V(:,:,s));
        level = 0.7 * max(I(:));
        
        L = I + level * BW;% * single(roi);

        imshow(L,[]);
        
        drawnow();
        
        %{
        if any(roi(:))
            pause();
            out_dir = fullfile('C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\qc-processing\temp', ts);
            if ~exist(out_dir, 'dir'), mkdir(out_dir); end
            saveas(hf, fullfile(out_dir, ['Image-' num2str(s, '%03d') '.png']));
        end
        %}

    end

end