


function util_make_null_mask_file()

    %BASE_DIR     = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\raw-data\7570-7862\Segmentation\7862';
    BASE_DIR     = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\raw-data\9310-9506\Segmentation\9442';

    NII_FILE_IN  = '9442_DCE_prone_PG.nii.gz';
    NII_FILE_OUT = '9442_DCE_prone_air.nii';

    H = niftiinfo(fullfile(BASE_DIR, NII_FILE_IN));
    V = niftiread(H);
    
    V = uint8(zeros(size(V)));
    
    niftiwrite(V, fullfile(BASE_DIR, NII_FILE_OUT), H, 'Compressed',true);

end